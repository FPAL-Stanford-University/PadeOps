module SolidGrid
    use kind_parameters,         only: rkind, clen
    use constants,               only: zero,eps,third,half,one,two,three,four,pi
    use FiltersMod,              only: filters
    use GridMod,                 only: grid
    use gridtools,               only: alloc_buffs, destroy_buffs
    use sgrid_hooks,             only: meshgen, initfields, hook_output, hook_bc, hook_timestep, hook_mixture_source
    use decomp_2d,               only: decomp_info, get_decomp_info, decomp_2d_init, decomp_2d_finalize, &
                                       transpose_x_to_y, transpose_y_to_x, transpose_y_to_z, transpose_z_to_y
    use DerivativesMod,          only: derivatives
    use DerivativesStaggeredMod, only: derivativesStagg
    use InterpolatorsMod,        only: interpolators
    use LADMod,                  only: ladobject
    use StiffGasEOS,             only: stiffgas
    use Sep1SolidEOS,            only: sep1solid
    use SolidMixtureMod,         only: solid_mixture
    use IOsgridMod,              only: IOsgrid

    use operators,               only: filter3D
   
    implicit none

    integer, parameter :: rho_index    = 1 
    integer, parameter :: u_index      = 2
    integer, parameter :: v_index      = 3
    integer, parameter :: w_index      = 4
    integer, parameter :: p_index      = 5
    integer, parameter :: T_index      = 6
    integer, parameter :: e_index      = 7
    integer, parameter :: sos_index    = 8
    integer, parameter :: mu_index     = 9
    integer, parameter :: bulk_index   = 10
    integer, parameter :: kap_index    = 11
    integer, parameter :: sxx_index    = 12
    integer, parameter :: sxy_index    = 13
    integer, parameter :: sxz_index    = 14
    integer, parameter :: syy_index    = 15
    integer, parameter :: syz_index    = 16
    integer, parameter :: szz_index    = 17
    integer, parameter :: rhoint_index = 18
    integer, parameter :: uint_index   = 19
    integer, parameter :: vint_index   = 20
    integer, parameter :: wint_index   = 21
    integer, parameter :: pint_index   = 22
    integer, parameter :: eint_index   = 23
    integer, parameter :: qyint_index  = 24
    integer, parameter :: tauxyint_index   = 25
    integer, parameter :: tauyyint_index   = 26
    integer, parameter :: tauyzint_index   = 27
    integer, parameter :: xflux_x_index       = 28
    integer, parameter :: xflux_y_index       = 29
    integer, parameter :: xflux_z_index       = 30
    integer, parameter :: yflux_x_index       = 31
    integer, parameter :: yflux_y_index       = 32
    integer, parameter :: yflux_z_index       = 33
    integer, parameter :: zflux_x_index       = 34
    integer, parameter :: zflux_y_index       = 35
    integer, parameter :: zflux_z_index       = 36
    integer, parameter :: xflux_e_index       = 37
    integer, parameter :: yflux_e_index       = 38
    integer, parameter :: zflux_e_index       = 39
    integer, parameter :: uJ_index            = 40
    integer, parameter :: vJ_index            = 41
    integer, parameter :: wJ_index            = 42
    integer, parameter :: keJ_index           = 43
    integer, parameter :: eJ_index            = 44
    integer, parameter :: tauxx_index         = 45
    integer, parameter :: tauyy_index         = 46
    integer, parameter :: tauzz_index         = 47
    integer, parameter :: tauxy_index         = 48
    integer, parameter :: tauyx_index         = 49
    integer, parameter :: tauxz_index         = 50
    integer, parameter :: tauzx_index         = 51
    integer, parameter :: tauyz_index         = 52
    integer, parameter :: tauzy_index         = 53
    integer, parameter :: tauxxe_index        = 54
    integer, parameter :: tauyye_index        = 55
    integer, parameter :: tauzze_index        = 56
    integer, parameter :: tauxye_index        = 57
    integer, parameter :: tauyxe_index        = 58
    integer, parameter :: tauyze_index        = 59
    integer, parameter :: tauzye_index        = 60
    integer, parameter :: tauxze_index        = 61
    integer, parameter :: tauzxe_index        = 62
    integer, parameter :: dudx_index          = 63
    integer, parameter :: dudy_index          = 64
    integer, parameter :: dudz_index          = 65
    integer, parameter :: dvdx_index          = 66
    integer, parameter :: dvdy_index          = 67
    integer, parameter :: dvdz_index          = 68
    integer, parameter :: dwdx_index          = 69
    integer, parameter :: dwdy_index          = 70
    integer, parameter :: dwdz_index          = 71 
    integer, parameter  :: tauSum_index       = 72
    integer, parameter  :: esum_index         = 73
    integer, parameter  :: esumJ_index        = 74
    integer, parameter  :: pmix_index         = 75
    integer, parameter  :: intP_index         = 76
    integer, parameter  :: qDiv_index         = 77
    integer, parameter  :: pEvolve_index      = 78
    integer, parameter  :: VFEvolve_index     = 79
    integer, parameter  :: pError_index       = 80
    integer, parameter  :: VFError_index      = 81
    integer, parameter  :: pJ_index           = 82
    integer, parameter  :: tauRho_index       = 83
    integer, parameter  :: metric_exact_index = 84
    integer, parameter  :: metric_index       = 85
    integer, parameter  :: metric_N2F_index   = 86
    integer, parameter  :: metric_half_index  = 87
    integer, parameter  :: nfields = 87

    integer, parameter :: mom_index = 1
    integer, parameter :: TE_index = mom_index+3
    integer, parameter :: ncnsrv  = TE_index

    ! These indices are for data management, do not change if you're not sure of what you're doing
    integer, parameter :: tauxyidx = 2
    integer, parameter :: tauxzidx = 3
    integer, parameter :: tauyzidx = 6
    integer, parameter :: tauxxidx = 4
    integer, parameter :: tauyyidx = 7
    integer, parameter :: tauzzidx = 8

    integer, parameter :: qxidx = 1
    integer, parameter :: qyidx = 5
    integer, parameter :: qzidx = 9


    ! Number of buffers to create
    integer, parameter :: nbufsx = 2
    integer, parameter :: nbufsy = 6
    integer, parameter :: nbufsz = 2


    type, extends(grid) :: sgrid
       
        type(filters),          allocatable :: gfil
        type(derivatives),      allocatable :: derD02, derD06,derD04,derCD06,der_nostretch,derCD06_nostretch
        type(solid_mixture),    allocatable :: mix
        type(ladobject),        allocatable :: LAD
        type(derivativesStagg), allocatable :: derStagg,derStagg_stretch
        type(derivativesStagg), allocatable :: derStaggd02
        type(interpolators),    allocatable ::interpMid
        type(interpolators),    allocatable ::interpMid02
        type( IOsgrid ),        allocatable :: viz

        logical     :: PTeqb                       ! Use pressure and temperature equilibrium formulation
        logical     :: pEqb                        ! nterpolators),    allocatable ::interpMidUse pressure equilibrium formulation
        logical     :: pRelax                      ! Use pressure and temperature non-equilibrium formulation, but relax pressure at each substep
        logical     :: use_gTg                     ! Use formulation with the Finger tensor g^T.g instead of the full g tensor
        logical     :: cnsrv_g, cnsrv_gt, cnsrv_gp, cnsrv_pe ! use conservative form of equations
        logical     :: strainHard                  ! use strainHardening
        logical     :: updateEtot                  ! Update species etot (vs ehydro) with pRelax
        logical     :: useOneG                     ! Use formulation with a single g or gTg field
        logical     :: useNC
        logical     :: intSharp                    ! Include interface sharpening terms
        logical     :: intSharp_cpl                ! Include coupling of sharpening with momentum and energy equations
        logical     :: intSharp_cpg                ! Include coupling of sharpening with kinematic equations
        logical     :: intSharp_cpg_west           ! Use form of kinematic sharpening terms derived by Jacob West
        logical     :: intSharp_spf                ! Use Shukla-Pantano-Freund method - not in divergence form
        logical     :: intSharp_ufv                ! Use finite volume discretization for sharpening term
        logical     :: intSharp_utw                ! Use Tiwari formulation
        logical     :: usePhiForm
        logical     :: twoPhaseLAD                 ! Use dYs/dx instead of the fickian 
        logical     :: LAD5eqn 
        real(rkind) :: intSharp_gam                ! Interface sharpening Gamma parameter
        real(rkind) :: intSharp_eps                ! Interface sharpening epsilon parameter
        real(rkind) :: intSharp_cut                ! Interface sharpening cutoff parameter, for VF approaching 1 or 0
        real(rkind) :: intSharp_dif                ! Interface sharpening VF out of bounds diffusion
        real(rkind) :: intSharp_tnh                ! Interface sharpening blending parameter
        real(rkind) :: intSharp_pfloor              ! Pressure floor for pressure-temperature relaxation / LAD
        real(rkind) :: intSharp_tfloor              ! Temperature floor for pressure-temperature relaxation / LAD
        logical :: intSharp_d02                    ! Use 2nd order
        logical :: intSharp_msk                    ! Mask FV diffusion
        logical :: intSharp_flt                    ! Use dealliasing filter for interface sharpening derivatives
        logical :: intSharp_flp                    ! Filter pressure

        logical     :: useAkshayForm	
        logical     :: weightedcurvature
	logical     :: surface_mask
        logical     :: LADInt,LADN2F
	logical     :: use_FV, use_D04, use_Stagg, use_XiLS         !flag to use FV in surface tension scheme
	logical     :: use_gradphi,energy_surfTen          !flag to use phi formulation in surface tension calculation
	logical     :: use_gradVF           !flag to use VF formulation in surface tension calculation		
        logical     :: use_gradXi
        logical     :: use_surfaceTension   !flag to turn on/off surface tension (in momentum and energy equations)
        logical     :: use_normFV
        logical     :: use_normInt
        logical     :: use_CnsrvSurfaceTension
        logical     :: Stretch1D, Stretch1Dy, Stretch1Dx, Stretch1Dz
        real(rkind) :: surfaceTension_coeff !constant coefficient for surface tension
        real(rkind) :: R, p_amb, XiLS_eps
        logical :: filt_mask = .FALSE.             ! mask filter in high gradient regions
        real(rkind), dimension(:,:,:,:), allocatable :: filt_tmp,filt_grad  ! temporary for filter mask
        real(rkind), dimension(:,:,:), allocatable :: filt_thrs  ! temporary for filter mask
        real(rkind) :: filt_cut  ! bulk threshold for filter mask

        real(rkind), dimension(:,:,:,:), allocatable :: Wcnsrv                               ! Conserved variables
        real(rkind), dimension(:,:,:,:), allocatable :: xbuf, ybuf, zbuf   ! Buffers
        real(rkind), dimension(:,:,:),   pointer     :: rho_int, u_int, v_int,w_int,tauxy_int, tauyy_int, tauyz_int,qy_int, e_int, TE, p_int
        real(rkind), dimension(:,:,:),   pointer     :: xflux_x, yflux_x, zflux_x, xflux_y, yflux_y, zflux_y, xflux_z, yflux_z, zflux_z, yflux_e, xflux_e, zflux_e

        real(rkind), dimension(:,:,:), pointer :: x 
        real(rkind), dimension(:,:,:), pointer :: y 
        real(rkind), dimension(:,:,:), pointer :: z 
       
        real(rkind), dimension(:,:,:), pointer :: eta1
        real(rkind), dimension(:,:,:), pointer :: eta2
        real(rkind), dimension(:,:,:), pointer :: eta3
 
        real(rkind), dimension(:,:,:), pointer :: rho 
        real(rkind), dimension(:,:,:), pointer :: u 
        real(rkind), dimension(:,:,:), pointer :: v 
        real(rkind), dimension(:,:,:), pointer :: w 
        real(rkind), dimension(:,:,:), pointer :: p 
        real(rkind), dimension(:,:,:), pointer :: T 
        real(rkind), dimension(:,:,:), pointer :: e 
        real(rkind), dimension(:,:,:), pointer :: sos 
        real(rkind), dimension(:,:,:), pointer :: mu 
        real(rkind), dimension(:,:,:), pointer :: bulk 
        real(rkind), dimension(:,:,:), pointer :: kap
        real(rkind), dimension(:,:,:), pointer :: tauaiidivu
        real(rkind), dimension(:,:,:), pointer :: pmix, intP
        real(rkind), dimension(:,:,:,:), pointer :: devstress
        real(rkind), dimension(:,:,:), pointer :: sxx
        real(rkind), dimension(:,:,:), pointer :: sxy
        real(rkind), dimension(:,:,:), pointer :: sxz
        real(rkind), dimension(:,:,:), pointer :: syy
        real(rkind), dimension(:,:,:), pointer :: syz
        real(rkind), dimension(:,:,:), pointer :: szz
       
        real(rkind), dimension(:,:,:), pointer :: tauxx, tauxxe, dudx
        real(rkind), dimension(:,:,:), pointer :: tauyy, tauyye, dudy
        real(rkind), dimension(:,:,:), pointer :: tauzz, tauzze, dudz
        real(rkind), dimension(:,:,:), pointer :: tauxy, tauxye, dvdx
        real(rkind), dimension(:,:,:), pointer :: tauyx, tauyxe, dvdy
        real(rkind), dimension(:,:,:), pointer :: tauxz, tauxze, dvdz
        real(rkind), dimension(:,:,:), pointer :: tauzx, tauzxe, dwdx
        real(rkind), dimension(:,:,:), pointer :: tauyz, tauyze,dwdy,metric_half, metric_N2F
        real(rkind), dimension(:,:,:), pointer :: tauzy, tauzye, dwdz,tauSum,esum, esumJ, metric, metric_exact
 
        real(rkind), dimension(:,:,:), pointer :: keJ, uJ, vJ, wJ, eJ, qDiv,pEvolve, VFEvolve, pError, VFerror, pJ, tauRho
        real(rkind) :: phys_mu1, phys_mu2
        real(rkind) :: phys_bulk1, phys_bulk2
        real(rkind) :: phys_kap1, phys_kap2
        real(rkind) :: st_limit, pthick,uthick,rhothick,Ys_wiggle,VF_wiggle,VF_thick,Ys_thick
        real(rkind), dimension(:,:,:,:), allocatable :: meshstretch
        real(rkind), dimension(:,:,:), allocatable :: yMetric,xMetric, zMetric,yMetric_half, xMetric_half, zMetric_half,yLADMetric, yMetric_F2N, dy_stretch
        contains
            procedure          :: init
            procedure          :: destroy
            procedure          :: laplacian
            procedure          :: gradient 
            procedure          :: secondder
            procedure          :: advance_RK45
            procedure          :: simulate
            procedure          :: update_p
            procedure          :: getRHS_P
            procedure          :: get_tauSum
            procedure          :: checkTau
            procedure, private :: get_dt
            procedure, private :: get_primitive
            procedure, private :: get_primitive_g
            procedure, private :: get_conserved
            procedure, private :: get_conserved_g
            procedure, private :: post_bc
            procedure, private :: post_bc_2
            procedure, private :: getRHS
            procedure, private :: getRHS_NC
            procedure, private :: getRHS_xStagg
            procedure, private :: getRHS_yStagg
            procedure, private :: getRHS_zStagg
            procedure, private :: getRHS_x
            procedure, private :: getRHS_y
            procedure, private :: getRHS_z
            procedure          :: filter
            procedure          :: getPhysicalProperties
            procedure, private :: get_tau
            procedure, private :: get_tauStagg
            procedure, private :: get_q
            procedure          :: get_qLAD
            procedure          :: coordinateTransform
    end type

contains
    subroutine init(this, inputfile )
        use reductions, only: P_MAXVAL
        use exits,      only: message, warning, nancheck, GracefulExit
        class(sgrid),target, intent(inout) :: this
        character(len=clen), intent(in) :: inputfile  

        integer :: nx, ny, nz
        integer :: ns = 1
        character(len=clen) :: outputdir
        character(len=clen) :: inputdir
        character(len=clen) :: vizprefix = "sgrid"
        real(rkind) :: tviz = zero
        character(len=clen), dimension(nfields) :: varnames
        logical :: periodicx = .true. 
        logical :: periodicy = .true. 
        logical :: periodicz = .true.
        character(len=clen) :: derivative_x = "cd10"  
        character(len=clen) :: derivative_y = "cd10" 
        character(len=clen) :: derivative_z = "cd10"
        character(len=clen) :: derivativeStagg_x = "d02" !default behavior for FV schemes is 2nd Order    
        character(len=clen) :: derivativeStagg_y = "d02" !default behavior for FV schemes is 2nd Order    
        character(len=clen) :: derivativeStagg_z = "d02" !default behavior for FV schemes is 2nd Order    
        character(len=clen) :: interpolator_x = "ei02"   !default behavior for FV schemes is 2nd Order    
        character(len=clen) :: interpolator_y = "ei02"   !default behavior for FV schemes is 2nd Order    
        character(len=clen) :: interpolator_z = "ei02"   !default behavior for FV schemes is 2nd Order    
        character(len=clen) :: filter_x = "cf90"  
        character(len=clen) :: filter_y = "cf90" 
        character(len=clen) :: filter_z = "cf90"
        integer :: prow = 0, pcol = 0 
        integer :: i, j, k, itmp(3) 
        integer :: ioUnit
        real(rkind) :: gam = 1.4_rkind
        real(rkind) :: q = 0
        real(rkind) :: Rgas = one
        real(rkind) :: PInf = zero
        real(rkind) :: shmod = zero
        integer :: nsteps = -1
        real(rkind) :: dt = -one
        real(rkind) :: tstop = one
        real(rkind) :: CFL = -one
        logical :: SkewSymm = .FALSE.
        real(rkind) :: Cmu = 0.002_rkind
        real(rkind) :: Cbeta = 0.5_rkind
        real(rkind) :: CbetaP = 0.0_rkind
        real(rkind) :: Ckap = 0.01_rkind
        real(rkind) :: CkapP = 0.0_rkind
        real(rkind) :: Cdiff = 0.003_rkind
        real(rkind) :: CY = 100._rkind
        real(rkind) :: Cvf1 = 0
        real(rkind) :: Cvf2 = 0
        real(rkind) :: Crho = 0
        real(rkind) :: Cdiff_g = 0.003_rkind
        real(rkind) :: Cdiff_gt = 0.003_rkind
        real(rkind) :: Cdiff_gp = 0.003_rkind
        real(rkind) :: Cdiff_pe = 0.003_rkind
        real(rkind) :: Cdiff_pe_2 = 100._rkind
        logical     :: useNC = .FALSE.
        logical     :: PTeqb = .TRUE., pEqb = .false., pRelax = .false., updateEtot = .false.
        logical     :: use_gTg = .FALSE., useOneG = .FALSE., intSharp = .FALSE., usePhiForm = .TRUE., intSharp_cpl = .TRUE., intSharp_cpg = .false., intSharp_cpg_west = .FALSE., intSharp_spf = .FALSE., intSharp_ufv = .TRUE., intSharp_utw = .FALSE., intSharp_d02 = .TRUE., intSharp_msk = .TRUE., intSharp_flt = .FALSE., intSharp_flp = .FALSE., strainHard = .FALSE., cnsrv_g = .FALSE., cnsrv_gt = .FALSE., cnsrv_gp = .FALSE., cnsrv_pe = .FALSE.
        logical     :: SOSmodel = .FALSE.      ! TRUE => equilibrium model; FALSE => frozen model, Details in Saurel et al. (2009)
        logical     :: useAkshayForm = .FALSE.,twoPhaseLAD = .FALSE.,LAD5eqn = .FALSE., use_CnsrvSurfaceTension = .FALSE., use_surfaceTension = .FALSE., use_normFV = .false., use_normInt = .false.,use_gradXi = .FALSE., use_gradphi = .FALSE., use_gradVF = .FALSE., use_Stagg = .FALSE., use_FV = .FALSE.,use_D04 = .FALSE., surface_mask = .FALSE., weightedcurvature = .FALSE. , energy_surfTen = .FALSE., use_XiLS = .FALSE., LADN2F = .FALSE., LADInt = .FALSE., Stretch1D = .FALSE., Stretch1Dx = .FALSE., Stretch1Dy = .FALSE., Stretch1Dz = .FALSE.
        real(rkind) :: surfaceTension_coeff = 0.0d0, R = 1d0, p_amb = 1d0 
        integer     :: x_bc1 = 0, x_bcn = 0, y_bc1 = 0, y_bcn = 0, z_bc1 = 0, z_bcn = 0    ! 0: general, 1: symmetric/anti-symmetric
        real(rkind) :: phys_mu1 = 0.0d0, phys_mu2 =0.0d0,Ys_wiggle = 0d0,VF_wiggle = 0d0, uthick = 0, rhothick = 0, pthick = 0,VF_thick = 0, Ys_thick = 0
        real(rkind) :: phys_bulk1 = 0.0d0, phys_bulk2 =0.0d0
        real(rkind) :: phys_kap1 = 0.0d0, phys_kap2 =0.0d0

        real(rkind) :: intSharp_gam = 0.0d0, intSharp_eps = 0.0d0, intSharp_cut = 1.0d-2, intSharp_dif = 1.0d1, intSharp_tnh = 1.0D-2, intSharp_pfloor = 0.0D0, intSharp_tfloor = 0.0D0, XiLS_eps = 0.0

        real(rkind) :: filter_alpha = 0.499

        namelist /INPUT/       nx, ny, nz, tstop, dt, CFL, nsteps, &
                             inputdir, outputdir, vizprefix, tviz, &
                                  periodicx, periodicy, periodicz, &
                         derivative_x, derivative_y, derivative_z, &
          derivativeStagg_x, derivativeStagg_y, derivativeStagg_z, &
                   interpolator_x, interpolator_y, interpolator_z, &
                                     filter_x, filter_y, filter_z, &
                                                       prow, pcol, &
                                                         SkewSymm, &                                                      
                                                     filter_alpha

        namelist /SINPUT/  gam, Rgas, PInf, shmod, &
                           useAkshayForm,twoPhaseLAD,LAD5eqn, useNC, PTeqb, pEqb, pRelax, SOSmodel, use_gTg, updateEtot, useOneG, intSharp, usePhiForm,intSharp_cpl, intSharp_cpg, intSharp_cpg_west, intSharp_spf, intSharp_ufv, intSharp_utw, intSharp_d02, intSharp_msk, intSharp_flt, intSharp_flp, intSharp_gam, intSharp_eps, intSharp_cut, intSharp_dif, intSharp_tnh, intSharp_pfloor, intSharp_tfloor, ns, Cmu, Cbeta, CbetaP, Ckap, CkapP,Cdiff, CY, Cvf1, Cvf2, Crho, Cdiff_g, Cdiff_gt, Cdiff_gp, Cdiff_pe, Cdiff_pe_2, &
                           x_bc1, x_bcn, y_bc1, y_bcn, z_bc1, z_bcn, &
                           strainHard, cnsrv_g, cnsrv_gt, cnsrv_gp, cnsrv_pe, phys_mu1, phys_mu2, phys_bulk1, phys_bulk2, phys_kap1, phys_kap2, &
                           use_CnsrvSurfaceTension, use_surfaceTension, use_normFV, use_normInt, use_gradXi,energy_surfTen,use_gradphi, use_gradVF, use_Stagg, use_FV,use_D04, surface_mask, weightedcurvature, &
                           surfaceTension_coeff, R, p_amb, use_XiLS,XiLS_eps,LADInt, LADN2F, Stretch1D, Stretch1Dx, Stretch1Dy, Stretch1Dz

        ioUnit = 11
        open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
        read(unit=ioUnit, NML=INPUT)
        read(unit=ioUnit, NML=SINPUT)
        close(ioUnit)


        this%nx = nx
        this%ny = ny
        this%nz = nz

        this%Ys_wiggle = Ys_wiggle
        this%VF_wiggle = VF_wiggle
        this%uthick    = uthick
        this%rhothick   = rhothick
        this%pthick     = pthick
        this%phys_mu1   = phys_mu1
        this%phys_mu2   = phys_mu2
        this%phys_bulk1 = phys_bulk1
        this%phys_bulk2 = phys_bulk2
        this%phys_kap1 = phys_kap1
        this%phys_kap2 = phys_kap2
        this%Ys_thick  = Ys_thick
        this%VF_thick  = VF_thick

        this%tsim = zero
        this%tstop = tstop
        this%dtfixed = dt
        this%dt = dt
        this%CFL = CFL

        this%step = 0
        this%nsteps = nsteps
        
        this%Stretch1Dy    = Stretch1Dy
        this%Stretch1Dx    = Stretch1Dx
        this%Stretch1Dz    = Stretch1Dz
        this%Stretch1D     = Stretch1D
        this%useAkshayForm = useAkshayForm
        this%twoPhaseLAD   = twoPhaseLAD
        this%LAD5eqn       = LAD5eqn
        this%PTeqb  = PTeqb
        this%pEqb   = pEqb
        this%pRelax = pRelax
        this%use_gTg = use_gTg
        this%updateEtot = updateEtot
        this%useOneG = useOneG
        this%useNC = useNC
        this%intSharp = intSharp
        this%intSharp_cpl = intSharp_cpl
        this%intSharp_cpg      = intSharp_cpg
        this%intSharp_cpg_west = intSharp_cpg_west
        this%intSharp_spf = intSharp_spf
        this%intSharp_ufv = intSharp_ufv
        this%intSharp_utw = intSharp_utw
        this%intSharp_gam = intSharp_gam
        this%intSharp_eps = intSharp_eps
        this%intSharp_cut = intSharp_cut
        this%intSharp_dif = intSharp_dif
        this%intSharp_tnh = intSharp_tnh
        this%intSharp_pfloor = intSharp_pfloor
        this%intSharp_tfloor = intSharp_tfloor
        this%intSharp_d02 = intSharp_d02
        this%intSharp_msk = intSharp_msk
        this%intSharp_flt = intSharp_flt
        this%intSharp_flp = intSharp_flp
        this%usePhiForm   = usePhiForm
        this%strainHard = strainHard
        this%cnsrv_g  = cnsrv_g
        this%cnsrv_gt = cnsrv_gt
        this%cnsrv_gp = cnsrv_gp
        this%cnsrv_pe = cnsrv_pe
         
        this%LADInt               = LADInt
        this%LADN2F               = LADN2F	
        this%weightedcurvature    = weightedcurvature
	this%surface_mask 	  = surface_mask
	this%use_FV  		  = use_FV
        this%use_Stagg            = use_Stagg
        this%use_D04              = use_D04
	this%use_gradphi 	  = use_gradphi
        this%energy_surfTen       = energy_surfTen
        this%use_gradXi           = use_gradXi
        this%use_normFV           = use_normFV
        this%use_normInt          = use_normInt
	this%use_gradVF 	  = use_gradVF
        this%use_surfaceTension   = use_surfaceTension 
        this%use_CnsrvSurfaceTension = use_CnsrvSurfaceTension 
        this%use_XiLS             = use_XiLS
        this%XiLS_eps             = XiLS_eps
        this%surfaceTension_coeff = surfaceTension_coeff
        this%R                    = R
        this%p_amb                = p_amb
        itmp(1:3) = 0; if(this%PTeqb) itmp(1) = 1; if(this%pEqb) itmp(2) = 1; if(this%pRelax) itmp(3) = 1; 

        if(sum(itmp) .ne. 1) then
            call GracefulExit("Exactly one among PTeqb, pEqb and pRelax should be true",4634)
        endif

        if(SOSmodel .and. this%pRelax) then
            call GracefulExit("Equilibrium sound speed model valid only with PTeqb or pEqb. SOSmodel must be .false. with pRelax",4634)
        endif

        if(this%useOneG .and. this%use_gTg) then
            call GracefulExit("Only g formulation supported with single g field for now", 4634)
        endif

        ! Allocate decomp
        if ( allocated(this%decomp) ) deallocate(this%decomp)
        allocate(this%decomp)
        
        ! Initialize decomp
        call decomp_2d_init(nx, ny, nz, prow, pcol)
        call get_decomp_info(this%decomp)
        ! Set all the attributes of the abstract grid type         
        this%outputdir = outputdir 
        
        this%periodicx = periodicx
        this%periodicy = periodicy
        this%periodicz = periodicz

        if ( ((x_bc1 /= 0) .AND. (x_bc1 /= 1)) ) call GracefulExit("x_bc1 can be only 0 (general) or 1 (symmetric)",4634)
        if ( ((x_bcn /= 0) .AND. (x_bcn /= 1)) ) call GracefulExit("x_bcn can be only 0 (general) or 1 (symmetric)",4634)
        if ( ((y_bc1 /= 0) .AND. (y_bc1 /= 1)) ) call GracefulExit("y_bc1 can be only 0 (general) or 1 (symmetric)",4634)
        if ( ((y_bcn /= 0) .AND. (y_bcn /= 1)) ) call GracefulExit("y_bcn can be only 0 (general) or 1 (symmetric)",4634)
        if ( ((z_bc1 /= 0) .AND. (z_bc1 /= 1)) ) call GracefulExit("z_bc1 can be only 0 (general) or 1 (symmetric)",4634)
        if ( ((z_bcn /= 0) .AND. (z_bcn /= 1)) ) call GracefulExit("z_bcn can be only 0 (general) or 1 (symmetric)",4634)
        this%x_bc = [x_bc1, x_bcn]
        this%y_bc = [y_bc1, y_bcn]
        this%z_bc = [z_bc1, z_bcn]

        this%derivative_x = derivative_x    
        this%derivative_y = derivative_y    
        this%derivative_z = derivative_z  

        this%filter_x = filter_x    
        this%filter_y = filter_y    
        this%filter_z = filter_z  
        
        ! Finally, set the local array dimensions
        this%nxp = this%decomp%ysz(1)
        this%nyp = this%decomp%ysz(2)
        this%nzp = this%decomp%ysz(3)

        ! Allocate mesh
        if ( allocated(this%mesh) ) deallocate(this%mesh) 
        call alloc_buffs(this%mesh,3,'y',this%decomp)

        ! Associate pointers for ease of use
        this%x    => this%mesh  (:,:,:, 1) 
        this%y    => this%mesh  (:,:,:, 2) 
        this%z    => this%mesh  (:,:,:, 3)

        ! Generate default mesh: X \in [-1, 1), Y \in [-1, 1), Z \in [-1, 1)
        this%dx = two/nx
        this%dy = two/ny
        this%dz = two/nz

        ! Generate default mesh 
        do k = 1,size(this%mesh,3)
            do j = 1,size(this%mesh,2)
                do i = 1,size(this%mesh,1)
                    this%mesh(i,j,k,1) = -one + (this%decomp%yst(1) - 1 + i - 1)*this%dx           
                    this%mesh(i,j,k,2) = -one + (this%decomp%yst(2) - 1 + j - 1)*this%dy           
                    this%mesh(i,j,k,3) = -one + (this%decomp%yst(3) - 1 + k - 1)*this%dz           
                end do 
            end do 
        end do  

        ! Go to hooks if a different mesh is desired 
        call meshgen(this%decomp, this%dx, this%dy, this%dz, this%mesh)

        if ( allocated(this%der_nostretch) ) deallocate(this%der_nostretch)
        allocate(this%der_nostretch)

         call this%der_nostretch%init(                           this%decomp, &
                              this%dx,       this%dy,        this%dz, &
                            periodicx,     periodicy,      periodicz, &
                            derivative_x,  derivative_y,   derivative_z, &
                              .false.,       .false.,        .false., &
                              .false.)

        if ( allocated(this%derCD06_nostretch) ) deallocate(this%derCD06_nostretch)
        allocate(this%derCD06_nostretch)

        call this%derCD06_nostretch%init(                           this%decomp, &
                             this%dx,       this%dy,        this%dz, &
                            periodicx,     periodicy,      periodicz, &
                            "cd06",        "cd06",         "cd06",    &
                              .false.,       .false.,        .false., &
                              .false.)

        print *, "dy = ", this%dy

        if ( allocated(this%derStagg_stretch) ) deallocate(this%derStagg_stretch)
        allocate(this%derStagg_stretch)

        call this%derStagg_stretch%init(                      this%decomp, &
                           this%dx,       this%dy,        this%dz, &
                         periodicx,     periodicy,      periodicz, &
          derivativeStagg_x, derivativeStagg_y, derivativeStagg_z, &
                           .false.,       .false.,        .false., &
                           .false.)

        if ( allocated(this%interpMid) ) deallocate(this%interpMid)
        allocate(this%interpMid)

        ! Initialize Interpolator
        call this%interpMid%init(                     this%decomp, &
                           this%dx,       this%dy,        this%dz, &
                         periodicx,     periodicy,      periodicz, &
                   interpolator_x, interpolator_y, interpolator_z, &
                           .false.,       .false.,        .false., &
                           .false.)
         if ( allocated(this%dy_stretch) ) deallocate(this%dy_stretch)
         allocate(this%dy_stretch(this%nxp,this%nyp,this%nzp) )
         if ( allocated(this%yMetric) ) deallocate(this%yMetric)
         allocate(this%yMetric(this%nxp,this%nyp,this%nzp) )
         if ( allocated(this%xMetric) ) deallocate(this%xMetric)
         allocate(this%xMetric(this%nxp,this%nyp,this%nzp) )
         if ( allocated(this%zMetric) ) deallocate(this%zMetric)
         allocate(this%zMetric(this%nxp,this%nyp,this%nzp) )
         if ( allocated(this%yLADMetric) ) deallocate(this%yLADMetric)
         allocate(this%yLADMetric(this%nxp,this%nyp,this%nzp) )
         if ( allocated(this%yMetric_F2N) ) deallocate(this%yMetric_F2N)
         allocate(this%yMetric_F2N(this%nxp,this%nyp,this%nzp) )

         if ( allocated(this%yMetric_half) ) deallocate(this%yMetric_half)
         allocate(this%yMetric_half(this%nxp,this%nyp,this%nzp) )
         if ( allocated(this%xMetric_half) ) deallocate(this%xMetric_half)
         allocate(this%xMetric_half(this%nxp,this%nyp,this%nzp) )
         if ( allocated(this%zMetric_half) ) deallocate(this%zMetric_half)
         allocate(this%zMetric_half(this%nxp,this%nyp,this%nzp) )

        if( this%Stretch1Dy ) then

            print *, "in Stretch 1D"
            if ( allocated(this%meshstretch) ) deallocate(this%meshstretch)
            call alloc_buffs(this%meshstretch,3,'y',this%decomp)
            print *, "allocation issue?"
            this%eta1    => this%meshstretch  (:,:,:, 1)
            this%eta2    => this%meshstretch  (:,:,:, 2)
            this%eta3  => this%meshstretch  (:,:,:, 3)

            call this%coordinateTransform(this%dx,this%dy,this%dz)
        endif
 
        ! Allocate der
        if ( allocated(this%der) ) deallocate(this%der)
        allocate(this%der)
        ! Allocate derD02
        if ( allocated(this%derD02) ) deallocate(this%derD02)
        allocate(this%derD02)
        
        if ( allocated(this%derD04) ) deallocate(this%derD04)
        allocate(this%derD04)
        ! Allocate derD06
        if ( allocated(this%derD06) ) deallocate(this%derD06)
        allocate(this%derD06)
         ! Allocate derCD06
        if ( allocated(this%derCD06) ) deallocate(this%derCD06)
        allocate(this%derCD06)
        ! Allocate derStagg
        if ( allocated(this%derStagg) ) deallocate(this%derStagg)
        allocate(this%derStagg)
        if ( allocated(this%derStaggd02) ) deallocate(this%derStaggd02)
        allocate(this%derStaggd02)

        if( this%Stretch1Dy) then
           call this%der%init(                           this%decomp, &
                           this%dx,       this%dy,        this%dz, &
                         periodicx,     periodicy,      periodicz, &
                      derivative_x,  derivative_y,   derivative_z, &
                           .false.,       .true.,        .false., &
                           .false.)

           call this%der%init_gridStretch1D(this%yMetric,this%yMetric_half,this%yLADMetric)

           ! Initialize derivatives
           call this%derD02%init(                           this%decomp, &
                              this%dx,       this%dy,        this%dz, &
                            periodicx,     periodicy,      periodicz, &
                                              "d02",  "d02",   "d02", &
                              .false.,       .true.,        .false., &
                              .false.)

          call this%derD02%init_gridStretch1D(this%yMetric,this%yMetric_half,this%yLADMetric)

          ! Initialize derivatives
           call this%derD04%init(                           this%decomp, &
                              this%dx,       this%dy,        this%dz, &
                            periodicx,     periodicy,      periodicz, &
                                              "d04",  "d04",   "d04", &
                              .false.,       .true.,        .false., &
                              .false.)

           call this%derD04%init_gridStretch1D(this%yMetric,this%yMetric_half,this%yLADMetric)

           ! Initialize derivatives
           call this%derD06%init(                           this%decomp, &
                              this%dx,       this%dy,        this%dz, &
                            periodicx,     periodicy,      periodicz, &
                                              "d06",  "d06",   "d06", &
                              .false.,       .true.,        .false., &
                              .false.)

           call this%derD06%init_gridStretch1D(this%yMetric,this%yMetric_half,this%yLADMetric)

           ! Initialize derivatives
           call this%derCD06%init(                           this%decomp, &
                              this%dx,       this%dy,        this%dz, &
                            periodicx,     periodicy,      periodicz, &
                                              "cd06",  "cd06",   "cd06", &
                              .false.,       .true.,        .false., &
                              .false.)
           call this%derCD06%init_gridStretch1D(this%yMetric,this%yMetric_half,this%yLADMetric)

                ! Initialize Staggered derivatives
           call this%derStagg%init(                      this%decomp, &
                           this%dx,       this%dy,        this%dz, &
                         periodicx,     periodicy,      periodicz, &
           derivativeStagg_x, derivativeStagg_y, derivativeStagg_z, &
                           .false.,       .true.,        .false., &
                           .false.)

           call this%derStagg%init_gridStretch1D(this%yMetric_F2N,this%yMetric_half,this%yLADMetric)

           ! Initialize Staggered derivatives
           call this%derStaggd02%init(                      this%decomp, &
                           this%dx,       this%dy,        this%dz, &
                         periodicx,     periodicy,      periodicz, &
                                              "d02", "d02", "d02", &
                           .false.,       .true.,        .false., &
                           .false.)

           call this%derStaggd02%init_gridStretch1D(this%yMetric_F2N,this%yMetric_half,this%yLADMetric)

        else 
        ! Initialize derivatives 
        call this%der%init(                           this%decomp, &
                           this%dx,       this%dy,        this%dz, &
                         periodicx,     periodicy,      periodicz, &
                      derivative_x,  derivative_y,   derivative_z, &
                           .false.,       .false.,        .false., &
                           .false.)      


        ! Initialize derivatives 
        call this%derD02%init(                           this%decomp, &
                              this%dx,       this%dy,        this%dz, &
                            periodicx,     periodicy,      periodicz, &
                                              "d02",  "d02",   "d02", &
                              .false.,       .false.,        .false., &
                              .false.)     



        
        call this%derD04%init(                           this%decomp, &
                              this%dx,       this%dy,        this%dz, &
                            periodicx,     periodicy,      periodicz, &
                                              "d04",  "d04",   "d04", &
                              .false.,       .false.,        .false., &
                              .false.)

 

 ! Initialize derivatives 
        call this%derD06%init(                           this%decomp, &
                              this%dx,       this%dy,        this%dz, &
                            periodicx,     periodicy,      periodicz, &
                                              "d06",  "d06",   "d06", &
                              .false.,       .false.,        .false., &
                              .false.)      

        ! Initialize derivatives 
        call this%derCD06%init(                           this%decomp, &
                              this%dx,       this%dy,        this%dz, &
                            periodicx,     periodicy,      periodicz, &
                                              "cd06",  "cd06",   "cd06", &
                              .false.,       .false.,        .false., &
                              .false.)


        ! Initialize Staggered derivatives 
        call this%derStagg%init(                      this%decomp, &
                           this%dx,       this%dy,        this%dz, &
                         periodicx,     periodicy,      periodicz, &
          derivativeStagg_x, derivativeStagg_y, derivativeStagg_z, &
                           .false.,       .false.,        .false., &
                           .false.)      



        ! Initialize Staggered derivatives
        call this%derStaggd02%init(                      this%decomp, &
                           this%dx,       this%dy,        this%dz, &
                         periodicx,     periodicy,      periodicz, &
                                              "d02", "d02", "d02", &
                           .false.,       .false.,        .false., &
                           .false.)

        end if

        print *, "allocated derivatives"
        ! Allocate interpMid
        if ( allocated(this%interpMid02) ) deallocate(this%interpMid02)
        allocate(this%interpMid02)

        ! Initialize Interpolator 
        call this%interpMid02%init(                     this%decomp, &
                           this%dx,       this%dy,        this%dz, &
                         periodicx,     periodicy,      periodicz, &
                                            "ei02", "ei02","ei02", &
                           .false.,       .false.,        .false., &
                           .false.)      

        


        ! Allocate fil and gfil
        if ( allocated(this%fil) ) deallocate(this%fil)
        allocate(this%fil)
        if ( allocated(this%gfil) ) deallocate(this%gfil)
        allocate(this%gfil)
        
        ! Initialize filters
        call this%fil%init(                           this%decomp, &
                         periodicx,     periodicy,      periodicz, &
                          filter_x,      filter_y,       filter_z ,filter_alpha )     
        call this%gfil%init(                          this%decomp, &
                         periodicx,     periodicy,      periodicz, &
                        "gaussian",    "gaussian",     "gaussian"  )      
        ! Allocate LAD object
        if ( allocated(this%LAD) ) deallocate(this%LAD)
        allocate(this%LAD)
        !call
        !this%LAD%init(this%decomp,this%der,this%gfil,2,this%dx,this%dy,this%dz,Cbeta,Cmu,Ckap,Cdiff,CY,Cdiff_g,Cdiff_gt,Cdiff_gp,Cdiff_pe,Cdiff_pe_2)
        call this%LAD%init(this%decomp,this%der,this%derStagg, this%interpMid,this%gfil,2,this%dx,this%dy,this%dz,Cbeta,CbetaP,Cmu,Ckap,CkapP,Cdiff,CY,Cdiff_g,Cdiff_gt,Cdiff_gp,Cdiff_pe,Cdiff_pe_2,Crho,Cvf1, Cvf2,Stretch1Dy)

        ! Allocate mixture
        if ( allocated(this%mix) ) deallocate(this%mix)
        allocate(this%mix)      
  


  ! Initialize derivatives
        call this%mix%init(this%decomp,this%der,this%derD02,this%derStagg,this%derStaggd02,this%derD06,this%derCD06,this%derD04,this%interpMid,this%interpMid02,this%use_Stagg,this%LADN2F,this%LADInt,this%fil,this%gfil,this%LAD,ns,this%PTeqb,this%pEqb,this%pRelax,SOSmodel,this%use_gTg,this%updateEtot,this%useAkshayForm,this%twoPhaseLAD,this%LAD5eqn,this%useOneG,this%intSharp,this%usePhiForm, this%intSharp_cpl,this%intSharp_cpg,this%intSharp_cpg_west,this%intSharp_spf,this%intSharp_ufv,this%intSharp_utw,this%intSharp_d02,this%intSharp_msk,this%intSharp_flt,this%intSharp_gam,this%intSharp_eps,this%intSharp_cut,this%intSharp_dif,this%intSharp_tnh,this%intSharp_pfloor,this%use_surfaceTension,this%use_normFV,this%use_normInt, this%use_gradXi,this%energy_surfTen, this%use_gradphi, this%use_gradVF, this%surfaceTension_coeff, this%use_FV,this%use_XiLS, this%XiLS_eps, this%use_D04, this%surface_mask, this%weightedcurvature, this%strainHard,this%cnsrv_g,this%cnsrv_gt,this%cnsrv_gp,this%cnsrv_pe,this%x_bc,this%y_bc,this%z_bc)

        !allocate(this%mix, source=solid_mixture(this%decomp,this%der,this%fil,this%LAD,ns))

        ! Allocate fields
        if ( allocated(this%fields) ) deallocate(this%fields) 
        call alloc_buffs(this%fields,nfields,'y',this%decomp)
        
        if ( allocated(this%Wcnsrv) ) deallocate(this%Wcnsrv) 
        call alloc_buffs(this%Wcnsrv,ncnsrv,'y',this%decomp)

        ! Allocate temporary filter variable
        if ( allocated(this%filt_tmp) ) deallocate(this%filt_tmp) 
        call alloc_buffs(this%filt_tmp,4+this%mix%ns,'y',this%decomp)

        ! Allocate temporary filter variable
        if ( allocated(this%filt_grad) ) deallocate(this%filt_grad) 
        call alloc_buffs(this%filt_grad,3,'y',this%decomp)

        ! Allocate temporary filter variable
        if ( allocated(this%filt_thrs) ) deallocate(this%filt_thrs) 
        allocate(this%filt_thrs(this%nxp,this%nyp,this%nzp) )


        ! Associate pointers for ease of use
        this%rho  => this%fields(:,:,:, rho_index) 
        this%u    => this%fields(:,:,:,   u_index) 
        this%v    => this%fields(:,:,:,   v_index) 
        this%w    => this%fields(:,:,:,   w_index)  
        this%p    => this%fields(:,:,:,   p_index)  
        this%T    => this%fields(:,:,:,   T_index)  
        this%e    => this%fields(:,:,:,   e_index)  
        this%sos  => this%fields(:,:,:, sos_index)  
        this%mu   => this%fields(:,:,:,  mu_index)  
        this%bulk => this%fields(:,:,:,bulk_index)  
        this%kap  => this%fields(:,:,:, kap_index)   
        this%pmix => this%fields(:,:,:,pmix_index) 
        this%devstress => this%fields(:,:,:,sxx_index:szz_index)
        this%sxx  => this%fields(:,:,:, sxx_index)   
        this%sxy  => this%fields(:,:,:, sxy_index)   
        this%sxz  => this%fields(:,:,:, sxz_index)   
        this%syy  => this%fields(:,:,:, syy_index)   
        this%syz  => this%fields(:,:,:, syz_index)   
        this%szz  => this%fields(:,:,:, szz_index)   
        this%rho_int => this%fields(:,:,:, rhoint_index)
        this%u_int => this%fields(:,:,:, uint_index)
        this%v_int => this%fields(:,:,:, vint_index)
        this%w_int => this%fields(:,:,:, wint_index)
        this%p_int => this%fields(:,:,:, pint_index)
        this%e_int => this%fields(:,:,:, eint_index)
        this%qy_int => this%fields(:,:,:, qyint_index)
        this%tauxy_int => this%fields(:,:,:, tauxyint_index)
        this%tauyy_int => this%fields(:,:,:, tauyyint_index)
        this%tauyz_int => this%fields(:,:,:, tauyzint_index)
        this%xflux_x   => this%fields(:,:,:, xflux_x_index)
        this%xflux_y   => this%fields(:,:,:, xflux_y_index)
        this%xflux_z   => this%fields(:,:,:, xflux_z_index)
        this%yflux_x   => this%fields(:,:,:, yflux_x_index)
        this%yflux_y   => this%fields(:,:,:, yflux_y_index)
        this%yflux_z   => this%fields(:,:,:, yflux_z_index)
        this%zflux_x   => this%fields(:,:,:, zflux_x_index)
        this%zflux_y   => this%fields(:,:,:, zflux_y_index)
        this%zflux_z   => this%fields(:,:,:, zflux_z_index)
        this%xflux_e   => this%fields(:,:,:, xflux_e_index)
        this%yflux_e   => this%fields(:,:,:, yflux_e_index)
        this%zflux_e   => this%fields(:,:,:, zflux_e_index)
        this%uJ        => this%fields(:,:,:, uJ_index)
        this%vJ        => this%fields(:,:,:, vJ_index)
        this%wJ        => this%fields(:,:,:, wJ_index)
        this%keJ       => this%fields(:,:,:, keJ_index)
        this%eJ        => this%fields(:,:,:,eJ_index)
        this%pJ        => this%fields(:,:,:,pJ_index)
        this%tauRho     => this%fields(:,:,:,tauRho_index) 
        this%tauxx     => this%fields(:,:,:,tauxx_index)
        this%tauyy     => this%fields(:,:,:,tauyy_index)
        this%tauzz     => this%fields(:,:,:,tauzz_index)
        this%tauxy     => this%fields(:,:,:,tauxy_index)
        this%tauyx     => this%fields(:,:,:,tauyx_index)
        this%tauxz     => this%fields(:,:,:,tauxz_index)
        this%tauzx     => this%fields(:,:,:,tauzx_index)
        this%tauyz     => this%fields(:,:,:,tauyz_index)
        this%tauzy     => this%fields(:,:,:,tauzy_index)
        this%tauxxe    => this%fields(:,:,:,tauxxe_index)
        this%tauyye    => this%fields(:,:,:,tauyye_index)
        this%tauzze    => this%fields(:,:,:,tauzze_index)
        this%tauxye    => this%fields(:,:,:,tauxye_index)
        this%tauyxe    => this%fields(:,:,:,tauyxe_index)
        this%tauyze    => this%fields(:,:,:,tauyze_index)
        this%tauzye    => this%fields(:,:,:,tauzye_index)
        this%tauxze    => this%fields(:,:,:,tauxze_index)
        this%tauzxe    => this%fields(:,:,:,tauzxe_index)
        this%dudx      => this%fields(:,:,:,dudx_index)
        this%dudy      => this%fields(:,:,:,dudy_index)
        this%dudz      => this%fields(:,:,:,dudz_index)
        this%dvdx      => this%fields(:,:,:,dvdx_index)
        this%dvdy      => this%fields(:,:,:,dvdy_index)
        this%dvdz      => this%fields(:,:,:,dvdz_index) 
        this%dwdx      => this%fields(:,:,:,dwdx_index)
        this%dwdy      => this%fields(:,:,:,dwdy_index) 
        this%dwdz      => this%fields(:,:,:,dwdz_index)
        this%tauSum    => this%fields(:,:,:,tauSum_index)
        this%esum      => this%fields(:,:,:,esum_index)
        this%esumJ     => this%fields(:,:,:,esumJ_index)
        this%intP      => this%fields(:,:,:,intP_index)
        this%qDiv      => this%fields(:,:,:,qDiv_index)
        this%pEvolve   => this%fields(:,:,:,pEvolve_index)
        this%VFEvolve  => this%fields(:,:,:,VFEvolve_index)
        this%pError    => this%fields(:,:,:,pError_index)
        this%VFError   => this%fields(:,:,:,VFError_index)       
        this%metric    => this%fields(:,:,:,metric_index)
        this%metric_exact    => this%fields(:,:,:,metric_exact_index)
        this%metric_half => this%fields(:,:,:,metric_half_index)
        this%metric_N2F  => this%fields(:,:,:,metric_N2F_index) 
        ! Initialize everything to a constant Zero
        this%fields = zero  

        ! Go to hooks if a different initialization is derired (Set mixture p, Ys, VF, u, v, w, rho)
        call initfields(this%decomp,this%der, this%dx, this%dy, this%dz, inputfile, this%mesh, this%fields, &
                        this%mix, this%tstop, this%dtfixed, tviz)
        ! Get hydrodynamic and elastic energies, stresses
        call this%mix%get_rhoYs_from_gVF(this%rho)  ! Get mixture rho and species Ys from species deformations and volume fractions
        call this%post_bc()
        ! Allocate 2 buffers for each of the three decompositions
        call alloc_buffs(this%xbuf,nbufsx,"x",this%decomp)
        call alloc_buffs(this%ybuf,nbufsy,"y",this%decomp)
        call alloc_buffs(this%zbuf,nbufsz,"z",this%decomp)

        this%SkewSymm = SkewSymm

        varnames( 1) = 'density'
        varnames( 2) = 'u'
        varnames( 3) = 'v'
        varnames( 4) = 'w'
        varnames( 5) = 'p'
        varnames( 6) = 'T'
        varnames( 7) = 'e'
        varnames( 8) = 'sos'
        varnames( 9) = 'mu'
        varnames(10) = 'bulk'
        varnames(11) = 'kap'
        varnames(12) = 'Sxx'
        varnames(13) = 'Sxy'
        varnames(14) = 'Sxz'
        varnames(15) = 'Syy'
        varnames(16) = 'Syz'
        varnames(17) = 'Szz'
        varnames(18) = 'rhoint'
        varnames(19) = 'uint'
        varnames(20) = 'vint'
        varnames(21) = 'wint'  
        varnames(22) = 'pint'
        varnames(23) = 'eint'
        varnames(24) = 'qyint'
        varnames(25) = 'tauxyint'
        varnames(26) = 'tauyyint' 
        varnames(27) = 'tauyzint'
        varnames(28) = 'xflux_x'
        varnames(29) = 'xflux_y'
        varnames(30) = 'xflux_z'
        varnames(31) = 'yflux_x'      
        varnames(32) = 'yflux_y' 
        varnames(33) = 'yflux_z'          
        varnames(34) = 'zflux_x'      
        varnames(35) = 'zflux_y' 
        varnames(36) = 'zflux_z'
        varnames(37) = 'xflux_e'      
        varnames(38) = 'yflux_e' 
        varnames(39) = 'zflux_e' 
        varnames(40) = 'uJ'
        varnames(41) = 'vJ'
        varnames(42) = 'wJ'
        varnames(43) = 'keJ'
        varnames(44) = 'eJ'  
        varnames(45) = 'tauxx'
        varnames(46) = 'tauyy'
        varnames(47) = 'tauzz'
        varnames(48) = 'tauxy'
        varnames(49) = 'tauyx'
        varnames(50) = 'tauxz'
        varnames(51) = 'tauzx'
        varnames(52) = 'tauyz'
        varnames(53) = 'tauzy'
        varnames(54) = 'tauxxe'
        varnames(55) = 'tauyye'
        varnames(56) = 'tauzze'
        varnames(57) = 'tauxye'
        varnames(58) = 'tauyxe'
        varnames(59) = 'tauyze'
        varnames(60) = 'tauzye'
        varnames(61) = 'tauxze'
        varnames(62) = 'tauzxe'
        varnames(63) = 'dudx'
        varnames(64) = 'dudy'
        varnames(65) = 'dudz'
        varnames(66) = 'dvdx'
        varnames(67) = 'dvdy'
        varnames(68) = 'dvdz'
        varnames(69) = 'dwdx'
        varnames(70) = 'dwdy'
        varnames(71) = 'dwdz'
        varnames(72) = 'tauSum'
        varnames(73) = 'esum'
        varnames(74) = "esumJ"
        varnames(75) = "pmix2"
        varnames(76) = "intP"
        varnames(77) = "qDiv"
        varnames(78) = "pEvolve"
        varnames(79) = "VFEvolve"
        varnames(80) = "pError"
        varnames(81) = "VFError"
        varnames(82) = "pJ"
        varnames(83) = 'tauRho'
        varnames(84) = 'metric_exact'
        varnames(85) = 'metric'
        varnames(86) = 'metric_N2F'
        varnames(87) = 'metric_half'
        allocate(this%viz)
        call this%viz%init(this%outputdir, vizprefix, nfields, varnames)
        this%tviz = tviz
    end subroutine


    subroutine destroy(this)
        class(sgrid), intent(inout) :: this

        ! Nullify pointers
        nullify(this%x); nullify(this%y); nullify(this%z)
        nullify(this%eta1); nullify(this%eta2); nullify(this%eta3)
        nullify(this%rho      )
        nullify(this%u        )
        nullify(this%v        )
        nullify(this%w        )
        nullify(this%p        )
        nullify(this%T        )
        nullify(this%e        )
        nullify(this%sos      )
        nullify(this%mu       )
        nullify(this%bulk     )
        nullify(this%kap      )
        nullify(this%devstress)
        nullify(this%sxx      )
        nullify(this%sxy      )
        nullify(this%sxz      )
        nullify(this%syy      )
        nullify(this%syz      )
        nullify(this%szz      )
        nullify(this%rho_int  )
        nullify(this%u_int    )
        nullify(this%v_int    )
        nullify(this%w_int    )
        nullify(this%p_int    )
        nullify(this%e_int    )
        nullify(this%qy_int   )
        nullify(this%tauxy_int)
        nullify(this%tauyy_int)
        nullify(this%tauyz_int)
        nullify(this%xflux_x  )
        nullify(this%xflux_y  )
        nullify(this%xflux_z  )
        nullify(this%yflux_x  )
        nullify(this%yflux_y  )
        nullify(this%yflux_z  )
        nullify(this%zflux_x  )
        nullify(this%zflux_y  ) 
        nullify(this%zflux_z  )
        nullify(this%xflux_e  )
        nullify(this%yflux_e  )
        nullify(this%zflux_e  )
        nullify(this%uJ)
        nullify(this%vJ)
        nullify(this%wJ)
        nullify(this%keJ)
        nullify(this%eJ)
        nullify(this%pJ)
        nullify(this%tauRho)
        nullify(this%tauxx)
        nullify(this%tauyy)
        nullify(this%tauzz)
        nullify(this%tauxy)
        nullify(this%tauyx)
        nullify(this%tauxz)
        nullify(this%tauzx)
        nullify(this%tauyz)
        nullify(this%tauzy)
        nullify(this%tauxxe)
        nullify(this%tauyye)
        nullify(this%tauzze)
        nullify(this%tauxye)
        nullify(this%tauyxe)
        nullify(this%tauyze)
        nullify(this%tauzye)
        nullify(this%tauxze)
        nullify(this%tauzxe)
        nullify(this%dudx)
        nullify(this%dudy)
        nullify(this%dudz)
        nullify(this%dvdx)
        nullify(this%dvdy)
        nullify(this%dvdz)
        nullify(this%dwdx)
        nullify(this%dwdy)
        nullify(this%dwdz)
        nullify(this%tauSum)
        nullify(this%esum)
        nullify(this%esumJ)
        nullify(this%pmix)
        nullify(this%intP)
        nullify(this%qDiv)
        nullify(this%pEvolve)
        nullify(this%VFEvolve)
        nullify(this%VFError)
        nullify(this%pError)
        nullify(this%metric)
        nullify(this%metric_exact)
        nullify(this%metric_half)
        nullify(this%metric_N2F)
        if (allocated(this%mesh)) deallocate(this%mesh) 
        if (allocated(this%fields)) deallocate(this%fields) 
        if (allocated(this%mesh)) deallocate(this%meshstretch)
        if (allocated(this%dy_stretch)) deallocate(this%dy_stretch)
        if (allocated(this%yMetric)) deallocate(this%yMetric)      
        if (allocated(this%yMetric_F2N)) deallocate(this%yMetric_F2N) 
        if (allocated(this%yLADMetric)) deallocate(this%yLADMetric) 
        if (allocated(this%xMetric)) deallocate(this%xMetric)        
        if (allocated(this%zMetric)) deallocate(this%zMetric)
        if (allocated(this%yMetric_half)) deallocate(this%yMetric_half)
        if (allocated(this%xMetric_half)) deallocate(this%xMetric_half)
        if (allocated(this%zMetric_half)) deallocate(this%zMetric_half)
        call this%der%destroy()
        if (allocated(this%der)) deallocate(this%der) 

        call this%derD02%destroy()
        if (allocated(this%derD02)) deallocate(this%derD02) 
    
        call this%der_nostretch%destroy()
        if (allocated(this%der_nostretch)) deallocate(this%der_nostretch)

        call this%derCD06_nostretch%destroy()
        if (allocated(this%derCD06_nostretch)) deallocate(this%derCD06_nostretch)
 
        call this%derD04%destroy()
        if (allocated(this%derD04)) deallocate(this%derD04)

        call this%derD06%destroy()
        if (allocated(this%derD06)) deallocate(this%derD06)
  
        call this%derCD06%destroy()
        if (allocated(this%derCD06)) deallocate(this%derCD06)

        call this%fil%destroy()
        if (allocated(this%fil)) deallocate(this%fil) 
        
        call this%gfil%destroy()
        if (allocated(this%gfil)) deallocate(this%gfil) 

        call this%derStagg%destroy()
        if (allocated(this%derStagg)) deallocate(this%derStagg) 

        call this%interpMid%destroy()
        if (allocated(this%interpMid)) deallocate(this%interpMid) 
       
        call this%derStaggd02%destroy()
        if (allocated(this%derStaggd02)) deallocate(this%derStaggd02)

        call this%interpMid02%destroy()
        if (allocated(this%interpMid02)) deallocate(this%interpMid02)

        call this%derStagg_stretch%destroy()
        if (allocated(this%derStagg_stretch)) deallocate(this%derStagg_stretch)
 
        call destroy_buffs(this%xbuf)
        call destroy_buffs(this%ybuf)
        call destroy_buffs(this%zbuf)

        if ( allocated(this%mix) ) deallocate(this%mix)

        if ( allocated(this%LAD) ) deallocate(this%LAD)
        
        if (allocated(this%Wcnsrv)) deallocate(this%Wcnsrv) 

        if (allocated(this%filt_tmp)) deallocate(this%filt_tmp) 

        if (allocated(this%filt_grad)) deallocate(this%filt_grad) 

        if (allocated(this%filt_thrs)) deallocate(this%filt_thrs) 
       
        call this%viz%destroy()
        if (allocated(this%viz)) deallocate(this%viz)

        call decomp_2d_finalize
        if (allocated(this%decomp)) deallocate(this%decomp) 

    end subroutine

    subroutine gradient(this, f, dfdx, dfdy, dfdz, x_bc, y_bc, z_bc)
        class(sgrid),target, intent(inout) :: this
        real(rkind), intent(in), dimension(this%nxp, this%nyp, this%nzp) :: f
        real(rkind), intent(out), dimension(this%nxp, this%nyp, this%nzp) :: dfdx
        real(rkind), intent(out), dimension(this%nxp, this%nyp, this%nzp) :: dfdy
        real(rkind), intent(out), dimension(this%nxp, this%nyp, this%nzp) :: dfdz
        integer, dimension(2), optional, intent(in) :: x_bc, y_bc, z_bc

        type(derivatives), pointer :: der
        type(decomp_info), pointer :: decomp
        real(rkind), dimension(:,:,:), pointer :: xtmp,xdum,ztmp,zdum
        
        der => this%der
        decomp => this%decomp
        xtmp => this%xbuf(:,:,:,1)
        xdum => this%xbuf(:,:,:,2)
        ztmp => this%zbuf(:,:,:,1)
        zdum => this%zbuf(:,:,:,2)

        ! Get Y derivatives
        call der%ddy(f,dfdy,y_bc(1),y_bc(2))

        ! Get X derivatives
        call transpose_y_to_x(f,xtmp,decomp)
        call der%ddx(xtmp,xdum,x_bc(1),x_bc(2))
        call transpose_x_to_y(xdum,dfdx)

        ! Get Z derivatives
        call transpose_y_to_z(f,ztmp,decomp)
        call der%ddz(ztmp,zdum,z_bc(1),z_bc(2))
        call transpose_z_to_y(zdum,dfdz)

    end subroutine 

    subroutine laplacian(this, f, lapf, x_bc, y_bc, z_bc)
        use timer
        class(sgrid),target, intent(inout) :: this
        real(rkind), intent(in), dimension(this%nxp, this%nyp, this%nzp) :: f
        real(rkind), intent(out), dimension(this%nxp, this%nyp, this%nzp) :: lapf
        integer, dimension(2), optional, intent(in) :: x_bc, y_bc, z_bc
        
        real(rkind), dimension(:,:,:), pointer :: xtmp,xdum,ztmp,zdum, ytmp
        type(derivatives), pointer :: der
        type(decomp_info), pointer :: decomp
        

        der => this%der
        decomp => this%decomp
        xtmp => this%xbuf(:,:,:,1)
        xdum => this%xbuf(:,:,:,2)
        ztmp => this%zbuf(:,:,:,1)
        zdum => this%zbuf(:,:,:,2)
        ytmp => this%ybuf(:,:,:,1)

        ! Get Y derivatives
        call der%d2dy2(f,lapf,y_bc(1),y_bc(2))
        
        ! Get X derivatives
        call transpose_y_to_x(f,xtmp,this%decomp) 
        call this%der%d2dx2(xtmp,xdum,x_bc(1),x_bc(2))
        call transpose_x_to_y(xdum,ytmp,this%decomp)

        lapf = lapf + ytmp

        ! Get Z derivatives
        call transpose_y_to_z(f,ztmp,this%decomp)
        call this%der%d2dz2(ztmp,zdum,z_bc(1),z_bc(2))
        call transpose_z_to_y(zdum,ytmp,this%decomp)
        
        lapf = lapf + ytmp

    end subroutine

    

    subroutine secondder(this, f, d2fdx2, d2fdy2, d2fdz2, x_bc, y_bc, z_bc)
        class(sgrid),target, intent(inout) :: this
        real(rkind), intent(in), dimension(this%nxp, this%nyp, this%nzp) :: f
        real(rkind), intent(out), dimension(this%nxp, this%nyp, this%nzp) :: d2fdx2
        real(rkind), intent(out), dimension(this%nxp, this%nyp, this%nzp) :: d2fdy2
        real(rkind), intent(out), dimension(this%nxp, this%nyp, this%nzp) :: d2fdz2
        integer, dimension(2), optional, intent(in) :: x_bc, y_bc, z_bc

        type(derivatives), pointer :: der
        type(decomp_info), pointer :: decomp
        real(rkind), dimension(:,:,:), pointer :: xtmp,xdum,ztmp,zdum

        der => this%der
        decomp => this%decomp
        xtmp => this%xbuf(:,:,:,1)
        xdum => this%xbuf(:,:,:,2)
        ztmp => this%zbuf(:,:,:,1)
        zdum => this%zbuf(:,:,:,2)

        ! Get Y derivative
        call this%der%d2dy2(f,d2fdy2,y_bc(1),y_bc(2))

        ! Get X derivative
        call transpose_y_to_x(f,xtmp,decomp)
        call this%der%d2dx2(xtmp,xdum,x_bc(1),x_bc(2))
        call transpose_x_to_y(xdum,d2fdx2,decomp)

        ! Get Z derivative
        call transpose_y_to_z(f,ztmp,decomp)
        call this%der%d2dz2(ztmp,zdum,z_bc(1),z_bc(2))
        call transpose_z_to_y(zdum,d2fdz2,decomp)

    end subroutine
    
    subroutine coordinateTransform(this,dx,dy,dz) 
       use constants,        only: one, half, pi
       use operators, only: divergence,gradient,divergenceFV,interpolateFV,interpolateFV_x,interpolateFV_y,interpolateFV_z,gradFV_x, gradFV_y, gradFV_z, gradFV_N2Fx, gradFV_N2Fy, gradFV_N2Fz
       class(sgrid), intent(inout) :: this
       real(rkind),  intent(inout) :: dx,dy,dz
       real(rkind), dimension(:,:,:), pointer :: x,y,z,eta1,eta2,eta3
       integer :: i,j,k
       integer :: nx, ny, nz, ix1, ixn, iy1, iyn, iz1, izn
       real(rkind) :: L, STRETCH_RATIO = 1.5, Lr, Lr_half
       real(rkind), dimension(this%nxp, this%nyp, this%nzp) :: y_half,eta2_half,tmpdy2,ymetric_half_exact
       real(rkind), dimension(this%nxp, this%nyp, this%nzp) :: eta2_int,tmp,tmp1,tmp2,tmp3, tmpeta, tmpeta2
       nx = this%decomp%xsz(1); ny = this%decomp%ysz(2); nz = this%decomp%zsz(3)

        ! If base decomposition is in Y
       ix1 = this%decomp%yst(1); iy1 = this%decomp%yst(2); iz1 = this%decomp%yst(3)
       ixn = this%decomp%yen(1); iyn = this%decomp%yen(2); izn = this%decomp%yen(3)

       L = 2.0 

       y_half = this%y + 0.5*this%dy
        
       this%eta1 = this%x
       this%eta3 = this%z
       this%xMetric = 1
       this%zMetric = 1
       this%xMetric_half = 1
       this%zMetric_half = 1
       tmpeta = atanh( 2.0*this%y / ( 1.0 + 1.0 / STRETCH_RATIO) )
       Lr =  L /( tmpeta(1, ny,1) - tmpeta(1,1,1))

       this%eta2 = Lr*tmpeta
       tmpeta2 = Lr*tmpeta
       call this%der_nostretch%ddy(tmpeta2, tmp2, this%y_bc(1),this%y_bc(2))
       this%yMetric = 1/tmp2
       call interpolateFV_y(this%decomp,this%interpMid,this%eta2,eta2_int,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        
       call gradFV_y(this%decomp,this%derStagg_stretch,eta2_int,tmp,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)     
       this%yMetric_F2N = 1 / tmp 
       call gradFV_N2Fy(this%decomp,this%derStagg_stretch,this%eta2,tmp3,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
       this%yMetric_half = 1.0 / tmp3
       this%yMetric_half(:,this%nyp,:) = 1.0
        
       call this%der_nostretch%d2dy2(this%eta2,tmpdy2,this%y_bc(1),this%y_bc(2))
       this%yLADMetric = 1.0 / tmpdy2
       this%dy_stretch(:,2:this%nyp-1,:) = abs((this%eta2(:,3:this%nyp,:)-this%eta2(:,1:this%nyp-2,:))/2)
       this%dy_stretch(:,1,:) = abs((this%eta2(:,2,:)-this%eta2(:,1,:)))
       this%dy_stretch(:,this%nyp,:) = abs((this%eta2(:,this%nyp,:)-this%eta2(:,this%nyp-1,:)) )
    end subroutine 

    subroutine update_P(this,Qtmpp,isub,dt,x,y,z,tsim,periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
        use decomp_2d,  only: nrank
        use RKCoeffs,   only: RK45_A,RK45_B
        class(sgrid), intent(inout) :: this
        integer, intent(in) :: isub
        real(rkind), intent(in) :: dt,tsim
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in)  :: x,y,z
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(inout)  :: Qtmpp
        integer, dimension(2), intent(in) :: x_bc, y_bc, z_bc
        logical :: periodicx,periodicy,periodicz
        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: rhsP ! for mass fraction equation
        real(rkind)  :: dx,dy

        dy = y(1,2,1) - y(1,1,1)

        call this%getRHS_P(rhsP,periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)

        ! advance sub-step
        if(isub==1) Qtmpp = zero                   ! not really needed since RK45_A(1) = 0
        Qtmpp = dt*rhsP + RK45_A(isub)*Qtmpp
        this%p =  this%p  + RK45_B(isub)*Qtmpp
    end subroutine

    subroutine getRHS_P(this,rhsP,periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
        use operators, only: divergence,gradient,divergenceFV,interpolateFV,interpolateFV_x,interpolateFV_y,interpolateFV_z, gradFV_x, gradFV_y, gradFV_z
        class(sgrid),                                         intent(inout)  :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(out)    :: rhsP
        integer, dimension(2), intent(in) :: x_bc, y_bc, z_bc
        logical :: periodicx,periodicy,periodicz
        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: tmp,tmp1,tmp2,tmp3,u_int,v_int, w_int,flux,divu, Gam,esum,ksum, esumJ, divup, tauSum, eta1, eta2, rhoSos1, rhoSos, EtaSum, rhoSos2
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,3) :: rho_int,Ys_int,rhoYs_int, p_int
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,2) :: sosm, rhom, eta
        integer :: i
       ! vcon = 100
        
         call this%get_tauSum(tauSum)
         rhsP = 0.0
         rhoSos = zero
         rhoSos2 = zero
         do i = 1,2
           call this%mix%material(i)%getLAD_VF(this%rho,periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
           call this%mix%material(i)%getSpeciesDensity(this%rho,rhom(:,:,:,i))
           call this%mix%material(i)%hydro%get_sos2(rhom(:,:,:,i),this%p,sosm(:,:,:,i))
           rhoSos = rhoSos + rhom(:,:,:,i)*sosm(:,:,:,i)/this%mix%material(i)%VF
           rhoSos2 = rhoSos2 + rhom(:,:,:,i)*sosm(:,:,:,i)*this%mix%material(i)%VF
         enddo

         eta(:,:,:,1) = (rhom(:,:,:,2)*sosm(:,:,:,2) - rhom(:,:,:,1)*sosm(:,:,:,1) ) / rhoSos
         eta(:,:,:,2) = (rhom(:,:,:,1)*sosm(:,:,:,1) - rhom(:,:,:,2)*sosm(:,:,:,2) ) / rhoSos

         Gam = 0
         esum = 0
         EtaSum = 0

         do i = 1,2

            Gam = Gam + this%mix%material(i)%VF/(this%mix%material(i)%hydro%gam- 1 )
            esum = esum + this%mix%material(i)%eh*this%mix%material(i)%rhom*(this%mix%material(i)%intSharp_aFV + this%mix%material(i)%intSharp_aDiffFV)
            EtaSum = EtaSum + this%mix%material(i)%eh*this%mix%material(i)%rhom*eta(:,:,:,i)
         enddo

         this%esum = esum

       if(.NOT. this%use_Stagg) then


         call this%gradient(this%p, tmp1, tmp2, tmp3, this%x_bc,  this%y_bc,this%z_bc)
         call divergence(this%decomp, this%der,this%u, this%v, this%w, divu, this%x_bc, this%y_bc, this%z_bc)
       
       else

         call interpolateFV_x(this%decomp,this%interpMid,this%u,u_int,periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
         call interpolateFV_y(this%decomp,this%interpMid,this%v,v_int,periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
         call interpolateFV_z(this%decomp,this%interpMid,this%w,w_int,periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
         call interpolateFV(this%decomp,this%interpMid,this%p,p_int,periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)

         !this%u_int = u_int
         !this%v_int = v_int
         !this%w_int = w_int
         !this%tauxy_int  = p_int(:,:,:,1)
         !this%tauyy_int = p_int(:,:,:,2)
         !this%tauyz_int = p_int(:,:,:,3)
         call gradFV_x(this%decomp,this%derStagg,-p_int(:,:,:,1),tmp1,periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
         call gradFV_y(this%decomp,this%derStagg,-p_int(:,:,:,2),tmp2,periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
         call gradFV_z(this%decomp,this%derStagg,-p_int(:,:,:,3),tmp3,periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
         call divergenceFV(this%decomp,this%derStagg,-u_int*p_int(:,:,:,1),-v_int*p_int(:,:,:,2),-w_int*p_int(:,:,:,3),divup,periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
         call divergenceFV(this%decomp,this%derStagg,-u_int,-v_int,-w_int,divu,periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)

       endif


         if( .NOT. this%intSharp) then
   
            !rhsP = tmp1*this%u + tmp2*this%v + tmp3*this%w + (1/Gam)*((this%rho*this%e +this%p)*divu +tauSum)
            rhsP = divup - this%p*divu + (1/Gam)*((this%rho*this%e +this%p)*divu + tauSum ) ! - this%tauRho) !+ EtaSum*divu) !+ this%rho*this%qDiv ) !+ EtaSum*divu )
            !rhsP = divup - this%p*divu + rhoSos2*divu + (1/Gam)*(tauSum + EtaSum*divu )
         else
            !rhsP = divup - this%p*divu +  (1/Gam)*((this%rho*this%e + this%p)*divu + this%mix%intSharp_pFV + tauSum) 
            !rhsP = tmp1*this%u + tmp2*this%v + tmp3*this%w + (1/Gam)*((this%rho*this%e +this%p)*divu + this%mix%intSharp_hFV - esum + tauSum)
            rhsP = divup - this%p*divu +  (1/Gam)*((this%rho*this%e +this%p)*divu + this%mix%intSharp_hFV + tauSum - esum ) ! - this%tauRho )! EtaSum*divu) !+ this%rho*this%qDiv) !+ EtaSum*divu )
            !rhsP = divup - this%p*divu +  rhoSos2*divu + (1/Gam)*( this%mix%intSharp_hFV - esum + tauSum + EtaSum*divu )
         endif


         if( this%twoPhaseLAD ) then

            esumJ = 0
            do i = 1,2
                
                call this%mix%material(i)%getLAD_VF(this%rho,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc) 
                esumJ = esumJ + this%mix%material(i)%eh*this%mix%material(i)%rhom*this%mix%material(i)%vfLAD 

            enddo

            !rhsP = rhsP + (1/Gam) * this%pJ 
            rhsP = rhsP + (1/Gam) * (this%eJ  - esumJ) 

         end if

         this%esumJ = esumJ
         this%intP = this%p

    end subroutine


    subroutine simulate(this)
        use reductions, only: P_MEAN
        use timer,      only: tic, toc
        use exits,      only: GracefulExit, message, check_exit
        use decomp_2d,  only: nrank
        class(sgrid), target, intent(inout) :: this

        logical :: tcond, vizcond, stepcond
        character(len=clen) :: stability
        real(rkind) :: cputime
        real(rkind), dimension(:,:,:,:), allocatable, target :: duidxj
        real(rkind), dimension(:,:,:), pointer :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
        real(rkind), dimension(:,:,:), pointer :: ehmix
        integer :: i, imat
      

        allocate( duidxj(this%nxp, this%nyp, this%nzp, 9) )
        ! Get artificial properties for initial conditions
        dudx => duidxj(:,:,:,1); dudy => duidxj(:,:,:,2); dudz => duidxj(:,:,:,3);
        dvdx => duidxj(:,:,:,4); dvdy => duidxj(:,:,:,5); dvdz => duidxj(:,:,:,6);
        dwdx => duidxj(:,:,:,7); dwdy => duidxj(:,:,:,8); dwdz => duidxj(:,:,:,9);
        
        call this%gradient(this%u, dudx, dudy, dudz, -this%x_bc,  this%y_bc,  this%z_bc)
        call this%gradient(this%v, dvdx, dvdy, dvdz,  this%x_bc, -this%y_bc,  this%z_bc)
        call this%gradient(this%w, dwdx, dwdy, dwdz,  this%x_bc,  this%y_bc, -this%z_bc)
        do i=1,this%mix%ns
            !if (this%use_gTg) then
                ! Project g tensor to SPD space
                !call this%mix%material(i)%elastic%make_tensor_SPD(this%mix%material(i)%g)
                !call this%mix%material(i)%elastic%make_tensor_SPD(this%mix%material(i)%g_t) !mca check
            !end if
            ! Get massfraction gradients in Ji
            call this%gradient(this%mix%material(i)%Ys,this%mix%material(i)%Ji(:,:,:,1),&
                               this%mix%material(i)%Ji(:,:,:,2),this%mix%material(i)%Ji(:,:,:,3), this%x_bc,  this%y_bc, this%z_bc)
        end do

        ! compute artificial shear and bulk viscosities
        call this%getPhysicalProperties()
        !call this%LAD%get_viscosities(this%rho,duidxj,this%mu,this%bulk,this%x_bc,this%y_bc,this%z_bc)
        call this%LAD%get_viscosities(this%rho,this%p,this%sos,duidxj,this%mu,this%bulk,this%x_bc,this%y_bc,this%z_bc,this%dt,this%intSharp_pfloor,this%yMetric,this%dy_stretch)
        if (this%PTeqb) then
            ehmix => duidxj(:,:,:,4) ! use some storage space
            ehmix = this%e
            do imat = 1, this%mix%ns
                ehmix = ehmix - this%mix%material(imat)%Ys * this%mix%material(imat)%eel
            enddo
            call this%LAD%get_conductivity(this%rho,this%p,ehmix,this%T,this%sos,this%kap,this%x_bc,this%y_bc,this%z_bc,this%intSharp_tfloor)
        end if

        ! compute species artificial conductivities and diffusivities
        call this%mix%getLAD(this%rho,this%p,this%e,this%u, this%v, this%w, duidxj,this%sos,this%yMetric,this%dy_stretch,this%use_gTg,this%strainHard,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc,this%intSharp_tfloor)  ! Compute species LAD (kap, diff, diff_g, diff_gt,diff_pe)
        nullify(dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,ehmix)
        deallocate( duidxj )

        ! ------------------------------------------------
        call this%get_dt(stability)
        !populate surface tension terms at initial condition
        if(this%use_surfaceTension) then
            if(this%mix%ns.ne.2) then
                call GracefulExit("Surface tension is not defined for single-species, and not implemented for more than 2 species",4634)
            endif

             call this%mix%get_surfaceTension(this%rho,this%x_bc,this%y_bc,this%z_bc,this%dx,this%dy,this%dz,this%periodicx,this%periodicy,this%periodicz,this%u,this%v,this%w)  ! Compute surface tension terms for momentum and energy equations
        !      call this%mix%get_gradp(this%rho,this%x_bc,this%y_bc,this%z_bc,this%dx,this%dy,this%dz,this%periodicx,this%periodicy,this%periodicz,this%u,this%v,this%w)

         !       this%mix%surfaceTension_f = this%mix%gradp
        endif

        !call this%mix%Test_Der_NP(this%x,this%y,this%z,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        !call this%mix%Test_1DStretch(this%x,this%y,this%z,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        this%metric = this%yMetric
        this%metric_exact = this%yLADMetric
        this%metric_N2F   = this%yMetric_F2N
        this%metric_half  = this%yMetric_half
        !call this%mix%Test_Der(this%x,this%y,this%z,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        !print *, "before test"
        !call this%mix%Test_Der_Periodic(this%x,this%y,this%z,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        if(this%use_CnsrvSurfaceTension) then
            if(this%mix%ns.ne.2) then
                call GracefulExit("Surface tension is not defined for single-species, and not implemented for more than 2 species",4634)
            endif

            call this%mix%get_surfaceTensionCnsrv(this%rho,this%x_bc,this%y_bc,this%z_bc,this%dx,this%dy,this%dz,this%periodicx,this%periodicy,this%periodicz,this%u,this%v,this%w)
            ! Compute surface tension terms for momentum and energy equations

        endif


       call this%mix%get_gradp(this%rho,this%x_bc,this%y_bc,this%z_bc,this%dx,this%dy,this%dz,this%periodicx,this%periodicy,this%periodicz,this%u,this%v,this%w)

       

        ! Write out initial conditions
        ! call hook_output(this%decomp, this%dx, this%dy, this%dz, this%outputdir, this%mesh, this%fields, this%mix, this%tsim, this%viz%vizcount)
        !call hook_output(this%decomp,this%der,this%dx,this%dy,this%dz,this%outputdir,this%mesh,this%fields,this%mix,this%tsim,this%viz%vizcount,this%pthick,this%uthick,this%rhothick,this%Ys_thick,this%VF_thick,this%Ys_wiggle,this%VF_wiggle,this%x_bc,this%y_bc,this%z_bc)

        if( this%Stretch1Dy) then
           print *, "pre vis"
           call this%viz%WriteViz(this%decomp, this%meshstretch, this%fields, this%mix, this%tsim)
           print *, "post vis"
        else
           call this%viz%WriteViz(this%decomp, this%mesh, this%fields, this%mix,this%tsim)
        endif

        vizcond = .FALSE.
        
        ! Check for visualization condition and adjust time step
        if ( (this%tviz > zero) .AND. (this%tsim + this%dt > this%tviz * this%viz%vizcount) ) then
            this%dt = this%tviz * this%viz%vizcount - this%tsim
            vizcond = .TRUE.
            stability = 'vizdump'
        end if

        tcond = .TRUE.
        ! Check tstop condition
        if ( (this%tstop > zero) .AND. (this%tsim >= this%tstop) ) then
            tcond = .FALSE.
        else if ( (this%tstop > zero) .AND. (this%tsim + this%dt >= this%tstop) ) then
            this%dt = this%tstop - this%tsim
        end if

        ! Check nsteps condition
        if ( (this%nsteps <= 0) .OR. (this%step < this%nsteps) ) then
            stepcond = .TRUE.
        else
            stepcond = .FALSE.
        end if

        if ( (this%tstop <= zero) .AND. (this%nsteps <= 0) ) then
            call GracefulExit('No stopping criterion set. Set either tstop or nsteps to be positive.', 345)
        end if

        ! Start the simulation while loop
        if(nrank==0) write(*,*) 'Starting time loop'
        do while ( tcond .AND. stepcond )
            ! Advance time
            call tic()
            call this%advance_RK45()
            call toc(cputime)

            !call this%mix%thick_calculations(this%rho, this%p,this%u,this%pthick,this%uthick,this%rhothick,this%Ys_thick,this%VF_thick,this%Ys_wiggle,this%VF_wiggle, this%dx)          
            call message(1,"Time",this%tsim)
            call message(1,"Step",this%step)
            call message(2,"Time step",this%dt)
            call message(2,"Stability limit: "//trim(stability))
            call message(2,"CPU time (in seconds)",cputime)
            call hook_timestep(this%decomp, this%mesh, this%fields, this%mix, this%step, this%tsim)
            ! Write out vizualization dump if vizcond is met 
           if (vizcond) then
                ! call hook_output(this%decomp, this%dx, this%dy, this%dz, this%outputdir, this%mesh, this%fields, this%mix, this%tsim, this%viz%vizcount)
                !call hook_output(this%decomp,this%der,this%dx,this%dy,this%dz,this%outputdir,this%mesh,this%fields,this%mix,this%tsim,this%viz%vizcount,this%pthick,this%uthick,this%rhothick,this%Ys_thick,this%VF_thick,this%Ys_wiggle,this%VF_wiggle,this%x_bc,this%y_bc,this%z_bc)

                if( this%Stretch1Dy) then               
                   call this%viz%WriteViz(this%decomp, this%meshstretch, this%fields, this%mix, this%tsim)
                else
                   call this%viz%WriteViz(this%decomp, this%mesh, this%fields,this%mix, this%tsim)
                endif
                vizcond = .FALSE.
           end if
            
            ! Get the new time step
            call this%get_dt(stability)
            ! Check for visualization condition and adjust time step
            if ( (this%tviz > zero) .AND. (this%tsim + this%dt >= this%tviz * this%viz%vizcount) ) then
                this%dt = this%tviz * this%viz%vizcount - this%tsim
                vizcond = .TRUE.
            end if

            ! Check tstop condition
            if ( (this%tstop > zero) .AND. (this%tsim >= this%tstop*(one - eps)) ) then
                tcond = .FALSE.
            else if ( (this%tstop > zero) .AND. (this%tsim + this%dt >= this%tstop*(one - eps)) ) then
                this%dt = this%tstop - this%tsim
                stability = 'stop'
                vizcond = .TRUE.
            end if

            ! Check nsteps condition
            if ( (this%nsteps <= 0) .OR. (this%step < this%nsteps) ) then
                stepcond = .TRUE.
            else
                stepcond = .FALSE.
            end if

            ! Check for exitpdo file
            if(check_exit(this%outputdir)) then
                ! call hook_output(this%decomp, this%dx, this%dy, this%dz, this%outputdir, this%mesh, this%fields, this%mix, this%tsim, this%viz%vizcount)
                !call hook_output(this%decomp,this%der,this%dx,this%dy,this%dz,this%outputdir,this%mesh,this%fields,this%mix,this%tsim,this%viz%vizcount,this%pthick,this%uthick,this%rhothick,this%Ys_thick,this%VF_thick,this%Ys_wiggle,this%VF_wiggle,this%x_bc,this%y_bc,this%z_bc)
                if( this%Stretch1Dy) then

                   call this%viz%WriteViz(this%decomp, this%meshstretch, this%fields, this%mix, this%tsim)

                else
                
                   call this%viz%WriteViz(this%decomp, this%mesh, this%fields,this%mix, this%tsim) 

                endif

                call GracefulExit("Found exitpdo file in working directory",1234)
            endif

        end do

    end subroutine

    subroutine advance_RK45(this)
        use RKCoeffs,   only: RK45_steps,RK45_A,RK45_B
        use exits,      only: message,nancheck,GracefulExit
        use reductions, only: P_MAXVAL, P_MINVAL
        use decomp_2d,  only: nrank
        use operators, only: divergence,gradient
        use constants,               only: pi
        class(sgrid), target, intent(inout) :: this

        real(rkind)                                               :: Qtmpt      ! Temporary variable for RK45
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,ncnsrv) :: rhs        ! RHS for conserved variables
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,ncnsrv) :: Qtmp             ! Temporary variable for RK45
        real(rkind), dimension(this%nxp,this%nyp,this%nzp)        :: divu,Qtmpp, pmix ! Velocity divergence for species energy eq
        real(rkind), dimension(this%nxp,this%nyp,this%nzp)        :: viscwork         ! Viscous work term for species energy eq
        real(rkind), dimension(this%nxp,this%nyp,this%nzp)        :: Fsource          ! Source term for possible use in VF, g eh eqns
        integer :: isub,i,j,k,l,imat,iter,ii,jj,kk
        real(rkind), dimension(:,:,:,:), allocatable, target :: duidxj
        real(rkind), dimension(:,:,:), pointer :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
        character(len=clen) :: charout


        allocate( duidxj(this%nxp, this%nyp, this%nzp, 9) )
        ! Get artificial properties for initial conditions
        dudx => duidxj(:,:,:,1); dudy => duidxj(:,:,:,2); dudz => duidxj(:,:,:,3);
        dvdx => duidxj(:,:,:,4); dvdy => duidxj(:,:,:,5); dvdz => duidxj(:,:,:,6);
        dwdx => duidxj(:,:,:,7); dwdy => duidxj(:,:,:,8); dwdz => duidxj(:,:,:,9);



        Qtmp  = zero
        Qtmpt = zero
        pmix  = zero
      

       do isub = 1,  RK45_steps


            if(this%use_CnsrvSurfaceTension) then
                call this%mix%get_surfaceTensionPE(this%rho,this%x_bc,this%y_bc,this%z_bc,this%dx,this%dy,this%dz,this%periodicx,this%periodicy,this%periodicz,this%u,this%v,this%w)
            endif

            call this%get_conserved()
            !!!!!!!!!!!!!!!!!!!!!! COMMENTED OUT G STUFF !!!!!!!!!!!!!!!!!!!!!!
            !call this%get_conserved_g()
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            if ( nancheck(this%Wcnsrv,i,j,k,l) ) then
                call message("Wcnsrv: ",this%Wcnsrv(i,j,k,l))
                !write(charout,'(A,I1,A,I5,A,4(I5,A))') "NaN encountered in solution (Wcnsrv) at substep ", isub, " of step ", this%step+1, " at (",i,", ",j,", ",k,", ",l,") of Wcnsrv"
                write(charout,'(A,I1,A,I5,A,4(I5,A))') "NaN encountered in solution (Wcnsrv) at substep ", isub, " of step ", this%step+1, " at (",i+this%decomp%yst(1)-1,", ",j+this%decomp%yst(2)-1,", ",k+this%decomp%yst(3)-1,", ",l,") of Wcnsrv"
                call GracefulExit(trim(charout), 999)
            end if
            call this%mix%checkNaN()
            ! Pre-compute stress, LAD, J, etc.
            ! call this%mix%getSOS(this%rho,this%p,this%sos)
            call this%gradient(this%u, dudx, dudy, dudz, -this%x_bc,this%y_bc,this%z_bc)
            call this%gradient(this%v, dvdx, dvdy, dvdz,  this%x_bc,-this%y_bc,this%z_bc)
            call this%gradient(this%w, dwdx, dwdy, dwdz,  this%x_bc,this%y_bc,-this%z_bc)
            call this%mix%getLAD(this%rho,this%p,this%e,this%u, this%v, this%w,duidxj,this%sos,this%yMetric,this%dy_stretch,this%use_gTg,this%strainHard,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc,this%intSharp_tfloor)  ! Compute species LAD (kap, diff, diff_g, diff_gt,diff_pe)
            call this%mix%get_J(this%rho)                                          ! Compute diffusive mass fluxes
            call this%mix%get_q(this%x_bc,this%y_bc,this%z_bc)                     ! Compute diffusive thermal fluxes (including enthalpy diffusion)
            if(this%intSharp) then
               if(this%mix%ns.ne.2) then
                  call GracefulExit("Problem if ns=1, should work for ns>2 but not tested",4634)
               endif

               ! if (this%step .LE. this%st_limit) then
               !    if ((isub.eq.one).and.(nrank.eq.0)) print *, "int clean"

               !endif
               do imat=1,this%mix%ns
                   this%mix%material(imat)%intSharp_a = zero
                   this%mix%material(imat)%intSharp_aDiff = zero
                   this%mix%material(imat)%intSharp_aFV = zero !-- problem?
                   this%mix%material(imat)%intSharp_R = zero
                   this%mix%material(imat)%intSharp_RDiff = zero
                   this%mix%material(imat)%intSharp_RFV = zero
               enddo
               this%mix%intSharp_f = zero
               this%mix%intSharp_fDiff = zero
               this%mix%intSharp_fFV = zero
               this%mix%intSharp_h = zero
               this%mix%intSharp_hDiff = zero
               this%mix%intSharp_hFV = zero
               this%mix%intSharp_kFV = zero

               call this%mix%get_intSharp(this%rho,this%x_bc,this%y_bc,this%z_bc,this%dx,this%dy,this%dz,this%periodicx,this%periodicy,this%periodicz,this%u,this%v,this%w)

            else      
                  ! !debug
                  ! if ((isub.eq.one).and.(nrank.eq.0)) print*,"overwriting intSharp"
                   do imat=1,this%mix%ns
                      this%mix%material(imat)%intSharp_a = zero
                      this%mix%material(imat)%intSharp_aDiff = zero
                      this%mix%material(imat)%intSharp_aFV = zero !-- problem?
                      this%mix%material(imat)%intSharp_R = zero
                      this%mix%material(imat)%intSharp_RDiff = zero
                      this%mix%material(imat)%intSharp_RFV = zero
                   enddo
                  this%mix%intSharp_f = zero
                  this%mix%intSharp_fDiff = zero
                  this%mix%intSharp_fFV = zero
                  this%mix%intSharp_h = zero
                  this%mix%intSharp_hDiff = zero
                  this%mix%intSharp_hFV = zero
                  this%mix%intSharp_kFV = zero
                  ! !end debug
                  

            endif

            !Calculate contributions of surface tension to all equation RHS's
            if(this%use_surfaceTension) then
                if(this%mix%ns.ne.2) then
                    call GracefulExit("Surface tension is not defined for single-species, and not implemented for more than 2 species",4634)
                endif

                call this%mix%get_surfaceTension(this%rho,this%x_bc,this%y_bc,this%z_bc,this%dx,this%dy,this%dz,this%periodicx,this%periodicy,this%periodicz,this%u,this%v,this%w)  ! Compute surface tension terms for momentum and energy equations

            endif

             if(this%use_CnsrvSurfaceTension) then
                if(this%mix%ns.ne.2) then
                    call GracefulExit("Surface tension is not defined forsingle-species, and not implemented for more than 2 species",4634)
                endif

                call this%mix%get_surfaceTensionCnsrv(this%rho,this%x_bc,this%y_bc,this%z_bc,this%dx,this%dy,this%dz,this%periodicx,this%periodicy,this%periodicz,this%u,this%v,this%w)
                ! Compute surface tension terms for momentum and energy equations

            endif

            ! Update total mixture conserved variables

            if (this%useNC) then
              call this%getRHS_NC(rhs,divu, viscwork)
            else
              call this%getRHS(rhs,divu,viscwork)
            endif
           !!!!!!!!!!!!!! UNCOMMENT            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           Qtmp  = this%dt*rhs  + RK45_A(isub)*Qtmp
           this%Wcnsrv = this%Wcnsrv + RK45_B(isub)*Qtmp
           !!!!!!!!!!!!!!!!!!! UNCOMMENT       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  

            ! calculate sources if they are needed
            if(.not. this%PTeqb) call this%mix%calculate_source(this%rho,divu,this%u,this%v,this%w,this%p,Fsource,this%x_bc,this%y_bc,this%z_bc) ! -- actually, source terms should be included for PTeqb as well --NSG

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! COMMENTED OUT UPDATE G             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Now update all the individual species variables
            !check eps
            !call this%mix%update_g(isub,max(this%dt,eps),this%rho,this%u,this%v,this%w,this%x,this%y,this%z,Fsource,this%tsim,this%x_bc,this%y_bc,this%z_bc)               ! g tensor

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            !call this%get_primitive_g()

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! UNCOMMENT             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            call this%mix%update_Ys(isub,this%dt,this%rho,this%u,this%v,this%w,this%x,this%y,this%z,this%tsim,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)               ! Volume Fraction
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! UNCOMENT             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !if (.NOT. this%PTeqb) then
            if(this%pEqb) then
                call this%mix%update_VF(isub,this%dt,this%rho,this%u,this%v,this%w,this%x,this%y,this%z,this%tsim,divu,Fsource,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)   
             !  call this%update_P(Qtmpp,isub,this%dt,this%x,this%y,this%z,this%tsim,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc) 
            elseif(this%pRelax) then
                call this%mix%update_VF(isub,this%dt,this%rho,this%u,this%v,this%w,this%x,this%y,this%z,this%tsim,divu,Fsource,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)                        ! Volume Fraction
                call this%mix%update_eh(isub,this%dt,this%rho,this%u,this%v,this%w,this%x,this%y,this%z,this%tsim,divu,viscwork,Fsource,this%devstress,this%x_bc,this%y_bc,this%z_bc) ! Hydrodynamic energy
            end if

            ! Integrate simulation time to keep it in sync with RK substep
            Qtmpt = this%dt + RK45_A(isub)*Qtmpt
            this%tsim = this%tsim + RK45_B(isub)*Qtmpt

!            this%p = 10
            !this%w = 0
            !this%e = 2.941
            !this%T = 10
            !this%u = 10
            !this%v = 0 !(sin(pi*this%y)**2)*sin(2*pi*this%x)*cos(pi*this%tsim/4)
            ! !print *, '-----', 8, this%Wcnsrv(179,1,1,1:4)
            ! !do i = 1, size(this%Wcnsrv,1)
            ! !    write(*,'(4(e21.14,1x))') this%Wcnsrv(i,1,1,1:4)
            ! !enddo

            ! this%filt_mask = .FALSE.!.TRUE.
            ! if(this%filt_mask) then
            !    if(isub.eq.RK45_steps) then
            !    if(this%step.gt.one) then
            !       if ((isub.eq.one).and.(nrank.eq.0)) print*,"masking filter"
            !       this%filt_tmp(:,:,:,1) = this%Wcnsrv(:,:,:,mom_index  )
            !       this%filt_tmp(:,:,:,2) = this%Wcnsrv(:,:,:,mom_index+1)
            !       this%filt_tmp(:,:,:,3) = this%Wcnsrv(:,:,:,mom_index+2)
            !       this%filt_tmp(:,:,:,4) = this%Wcnsrv(:,:,:,TE_index   )
            !       do imat=1,this%mix%ns
            !          this%filt_tmp(:,:,:,4+imat) = this%mix%material(imat)%consrv(:,:,:,1)
            !       enddo
                  
            !       call this%filter(this%filt_tmp(:,:,:,1), this%fil, 1,-this%x_bc, this%y_bc, this%z_bc)
            !       call this%filter(this%filt_tmp(:,:,:,2), this%fil, 1, this%x_bc,-this%y_bc, this%z_bc)
            !       call this%filter(this%filt_tmp(:,:,:,3), this%fil, 1, this%x_bc, this%y_bc,-this%z_bc)
            !       call this%filter(this%filt_tmp(:,:,:,4), this%fil, 1, this%x_bc, this%y_bc, this%z_bc)
                  
            !       ! Filter the individual species variables
            !       do imat=1,this%mix%ns
            !          call filter3D(this%decomp, this%fil, this%filt_tmp(:,:,:,4+imat),1,this%x_bc,this%y_bc,this%z_bc)
            !       enddo
                  
            !       this%filt_cut = 1.0D10!1.0D-2
                  
            !       this%filt_thrs = zero
            !       call this%gradient(this%rho,this%filt_grad(:,:,:,1),this%filt_grad(:,:,:,2),this%filt_grad(:,:,:,3), this%x_bc,  this%y_bc,  this%z_bc)
            !       this%filt_thrs = max(this%filt_thrs,sqrt(this%filt_grad(:,:,:,1)**two + this%filt_grad(:,:,:,2)**two + this%filt_grad(:,:,:,3)**two))
            !       call this%gradient(this%p,this%filt_grad(:,:,:,1),this%filt_grad(:,:,:,2),this%filt_grad(:,:,:,3), this%x_bc,  this%y_bc,  this%z_bc)
            !       this%filt_thrs = max(this%filt_thrs,sqrt(this%filt_grad(:,:,:,1)**two + this%filt_grad(:,:,:,2)**two + this%filt_grad(:,:,:,3)**two))
            !       call this%gradient(this%T,this%filt_grad(:,:,:,1),this%filt_grad(:,:,:,2),this%filt_grad(:,:,:,3), this%x_bc,  this%y_bc,  this%z_bc)
            !       this%filt_thrs = max(this%filt_thrs,sqrt(this%filt_grad(:,:,:,1)**two + this%filt_grad(:,:,:,2)**two + this%filt_grad(:,:,:,3)**two))
            !       call this%gradient(sqrt(this%u**two+this%v**two+this%w**two),this%filt_grad(:,:,:,1),this%filt_grad(:,:,:,2),this%filt_grad(:,:,:,3), this%x_bc,  this%y_bc,  this%z_bc)
            !       this%filt_thrs = max(this%filt_thrs,sqrt(this%filt_grad(:,:,:,1)**two + this%filt_grad(:,:,:,2)**two + this%filt_grad(:,:,:,3)**two))

            !       where(this%filt_thrs.lt.this%filt_cut) ! only update where low bulk
            !          this%Wcnsrv(:,:,:,mom_index  ) = this%filt_tmp(:,:,:,1)
            !          this%Wcnsrv(:,:,:,mom_index+1) = this%filt_tmp(:,:,:,2)
            !          this%Wcnsrv(:,:,:,mom_index+2) = this%filt_tmp(:,:,:,3)
            !          this%Wcnsrv(:,:,:,TE_index   ) = this%filt_tmp(:,:,:,4)
            !       endwhere
            !       do imat=1,this%mix%ns
            !h          where(this%filt_thrs.lt.this%filt_cut) ! only update where low bulk
            !             this%mix%material(imat)%consrv(:,:,:,1) = this%filt_tmp(:,:,:,4+imat)
            !          endwhere
            !       enddo
               
            !    endif
            !    endif
            ! else
            

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! UNCOMMENT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!              if(.NOT. this%use_Stagg) then
                  ! Filter the conserved variables
                  !call this%filter(this%Wcnsrv(:,:,:,mom_index  ), this%fil, 1,-this%x_bc, this%y_bc, this%z_bc)
                  !call this%filter(this%Wcnsrv(:,:,:,mom_index+1), this%fil, 1, this%x_bc,-this%y_bc, this%z_bc)
                  !call this%filter(this%Wcnsrv(:,:,:,mom_index+2), this%fil, 1, this%x_bc, this%y_bc,-this%z_bc)
                  !call this%filter(this%Wcnsrv(:,:,:, TE_index  ), this%fil, 1, this%x_bc, this%y_bc, this%z_bc)
                  !call this%filter(this%p, this%fil, 1,this%x_bc, this%y_bc, this%z_bc)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! UNCOMMENT  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                  !do i = 1, size(this%Wcnsrv,1)
                  !    write(*,'(4(e21.14,1x))') this%Wcnsrv(i,1,1,1:4)
                  !enddo
             
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! UNCOMMENT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  ! Filter the individual species variables
                  !call this%mix%filter(1, this%x_bc, this%y_bc, this%z_bc)
!              end if 
           !  endif
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! UNCOMMENT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            if(this%use_CnsrvSurfaceTension) then

               call this%mix%getYs(this%rho)
               call this%mix%get_surfaceTensionPE(this%rho,this%x_bc,this%y_bc,this%z_bc,this%dx,this%dy,this%dz,this%periodicx,this%periodicy,this%periodicz,this%u,this%v,this%w)

            endif

            call this%get_primitive()
            !!!!!!!!!!! COMMENTED OUT THINGS TO DO WITH G             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !call this%get_primitive_g()
            !call this%mix%implicit_plastic(this%rho) !implicit plastic deformation using new g and new rho
            !call this%mix%filter_g(1, this%x_bc, this%y_bc, this%z_bc)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! if (.NOT. this%explPlast) then
            !     if (this%plastic) then
            !         ! Effect plastic deformations
            !         ! call this%elastic%plastic_deformation(this%g)
            !         call this%mix%plastic_deformation()
            !         call this%get_primitive()

            !         ! Filter the conserved variables
            !         do i = 1,5
            !             call this%filter(this%Wcnsrv(:,:,:,i), this%fil, 1)
            !         end do
            !         ! Filter the g tensor
            !         do i = 1,9
            !             call this%filter(this%g(:,:,:,i), this%fil, 1)
            !         end do
            !     end if
            ! end if
            
            if (this%PTeqb) then
               !call this%mix%equilibratePressureTemperature(this%rho, this%e, this%p, this%T, isub)
               ! do i=1,2

               !    call this%mix%equilibrateTemperature(this%rho, this%e, this%p, this%T, isub, RK45_steps) 
               !  call this%mix%equilibratePressureTemperature_new(this%rho, this%e, this%p, this%T, isub, RK45_steps) !fixes problem when negative mass fraction
               !     call this%mix%pressureLiquidGas(this%rho, this%e, this%p)
               !    call this%mix%get_pmix(this%p)                         ! Get mixture pressure
               !    call this%mix%get_Tmix(this%T)                         ! Get mixture temperature
               ! enddo

            elseif (this%pEqb) then
               call this%mix%equilibratePressure(this%rho, this%e, pmix)
               !call this%mix%updateP_VF(this%rho,this%e,this%p)
            elseif (this%pRelax) then
                call this%mix%relaxPressure(this%rho, this%e, this%p)
                !call this%mix%relaxPressure_os(this%rho, this%u, this%v, this%w, this%e, this%dt, this%p)
            end if
            
            !this%pError = abs(this%pEvolve - this%p)
            !this%VFError = abs(this%VFEvolve - this%mix%material(1)%VF)

            !print *, nrank, 11
            !#####################################################################
            if(this%intSharp_flp) then  !high pressure ratio shocks require additional dealiasing
               call this%filter(this%p, this%fil, 1, this%x_bc, this%y_bc, this%z_bc)
               call this%filter(this%T, this%fil, 1, this%x_bc, this%y_bc, this%z_bc)
               ! call this%filter(this%u, this%fil, 1,-this%x_bc, this%y_bc, this%z_bc)
               ! call this%filter(this%v, this%fil, 1, this%x_bc,-this%y_bc, this%z_bc)
               ! call this%filter(this%w, this%fil, 1, this%x_bc, this%y_bc,-this%z_bc)
            endif
            !#####################################################################

            
            call hook_bc(this%decomp, this%mesh, this%fields, this%mix, this%tsim, this%x_bc, this%y_bc, this%z_bc)
            if(this%pEqb) then
            call this%post_bc_2()
            else
            call this%post_bc()
            endif
            !call hook_output(this%decomp,this%der,this%dx,this%dy,this%dz,this%outputdir,this%mesh,this%fields,this%mix,this%tsim,this%viz%vizcount,this%x_bc,this%y_bc,this%z_bc)      
        end do

          
      !  if(this%pEqb) then
      !    call this%mix%MLong_relax(this%rho, this%e, this%p)
      !    call this%post_bc_2()
      !  endif

        this%step = this%step + 1
        nullify(dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz)
        deallocate( duidxj )
 
    end subroutine

    subroutine get_dt(this,stability)
        use reductions, only : P_MAXVAL, P_MINVAL
        use decomp_2d,  only: nrank
	use constants,        only: zero,third,half,twothird,one,two,seven,pi,eps
        class(sgrid), target, intent(inout) :: this
        character(len=*), intent(out) :: stability
        real(rkind) :: dtCFL, a, dtsigma,dtsigma2,dtsigma3,dtmu, dtbulk, dtkap, dtdiff, dtdiff_g, dtdiff_gt, dtdiff_gp, dtplast, phys_mu, delta, dtSharp_diff, dtSharp_Adiff,alpha,dtSharp_bound,st_fac=10.D0,dtYs1, dtYs2,dtVF1,dtVF2,deltay
        integer :: i
        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: uk   !Source term for possible use in VF, g eh eqns
        character(len=30) :: str,str2

        this%st_limit = 20

        if(this%Stretch1Dy) then
          deltay = P_MINVAL(this%dy_stretch)
          delta = min(this%dx,deltay,this%dz)
        else 
          delta = min(this%dx, this%dy, this%dz)
	endif
	phys_mu = min(this%phys_mu1, this%phys_mu2)

        ! continuum
        !dtCFL  = this%CFL / P_MAXVAL( ABS(this%u)/this%dx + ABS(this%v)/this%dy + ABS(this%w)/this%dz   &
        !       + this%sos*sqrt( one/(this%dx**2) + one/(this%dy**2) + one/(this%dz**2) ))

        dtCFL  = this%CFL / P_MAXVAL( ABS(this%u)/this%dx + ABS(this%v)/this%dy +  &
               + this%sos*sqrt( one/(this%dx**2) + one/(this%dy**2)  ))
        !print *, "u maxval"
        !print *, P_MAXVAL( ABS(this%u)/this%dx )
        !print *, "v maxval"
        !print *, P_MAXVAL( ABS(this%v)/this%dy )
        !print *, "w maxval"
        !print *, P_MAXVAL( ABS(this%w)/this%dz )        
        !print *, "sos"
        !print *, P_MAXVAL( this%sos*sqrt( one/(this%dx**2) + one/(this%dy**2) +one/(this%dz**2) ) )
        !print *, "dy"
        !print *, this%dy
        !print *, "dx"
        !print *, this%dx
        !print *, "dz"
        !print *, this%dz
        !print *, "sosmax"
        !print *, P_MAXVAL( this%sos)
        dtmu   = 0.2_rkind * delta**2 / (P_MAXVAL( this%mu  / this%rho ) + eps)! * this%CFL
        dtYs1 = 0.75_rkind * delta**2 / (P_MAXVAL( this%mix%material(1)%rhodiff  / this%rho ) + eps) 
        dtYs2 = 0.75_rkind * delta**2 / (P_MAXVAL( this%mix%material(2)%rhodiff / this%rho ) + eps)
        dtVF1 = 0.75_rkind * delta**2 / (P_MAXVAL( this%mix%material(1)%adiff / this%rho ) + eps)
        dtVF2 = 0.75_rkind * delta**2 / (P_MAXVAL( this%mix%material(2)%adiff /this%rho ) + eps)

        !dtbulk = 0.2_rkind * delta**2 / (P_MAXVAL( this%bulk/ this%rho ) + eps) * this%CFL
        dtbulk = 0.2_rkind * delta**2 / (P_MAXVAL( this%bulk/ this%rho ) + eps) !/ 5.0 !test /5
	
	if ((this%use_surfaceTension) .OR. (this%use_CnsrvSurfaceTension)) then
              !  if ( phys_mu > eps) then
          ! filter3D(this%decomp,this%fil,gradphi(:,:,:,3),iflag,x_bc,y_bc,z_bc)  
                  uk = ABS(this%u*this%mix%norm(:,:,:,1)) + ABS(this%v*this%mix%norm(:,:,:,2)) + ABS(this%w*this%mix%norm(:,:,:,3))
	          dtsigma = 4*max(1/P_MAXVAL(sqrt((4*pi*this%surfaceTension_coeff*(1/this%dx**3 + 1/this%dy**3 ))/ (this%rho))) , &
                                  P_MINVAL(this%mu/(this%surfaceTension_coeff*(1/this%dx + 1/this%dy ) ) ))
                                  !max(1/P_MAXVAL( sqrt((4*pi*this%surfaceTension_coeff*(1/this%dx**3 + 1/this%dy**3 +1/this%dx**3) )/ (this%rho))) , &
                                  !P_MINVAL(this%mu/(this%surfaceTension_coeff*(1/this%dx + 1/this%dy + 1/this%dz) ) ))
              !  else 
                  dtsigma2 = P_MINVAL( sqrt( (this%rho)/(4*pi*this%surfaceTension_coeff*(1/this%dx**3 + 1/this%dy**3 + 1/this%dx**3) ) ) )
                  dtsigma3 =  P_MINVAL( 1/ (sqrt((4*pi*this%surfaceTension_coeff) / (this%rho*( this%dx**3 + this%dy**3 +this%dz**3))) + uk*(1/this%dx + 1/this%dy + 1/this%dz  )))

              !  end if
	end if
        ! species specific
        call this%mix%get_dt(this%rho, delta, dtkap, dtdiff, dtdiff_g, dtdiff_gt, dtdiff_gp, dtplast)

        if (this%PTeqb) then
            !dtkap  = delta**2 / (P_MAXVAL( this%kap*this%T/(this%rho*this%sos**2)) + eps)   ! Cook (2007) formulation
            dtkap  = one / ( (P_MAXVAL(this%kap*this%T/(this%rho*delta**4)))**(third) + eps) ! Cook (2009) formulation
        end if

        dtkap     = 0.2_rkind * dtkap! * this%CFL
        !dtkap     = 0.2_rkind * dtkap !/ 5.0! test

        dtdiff    = 0.2_rkind * dtdiff! * this%CFL
        dtdiff_g  = 0.2_rkind * dtdiff_g! * this%CFL
        dtdiff_gt = 0.2_rkind * dtdiff_gt! * this%CFL
        dtdiff_gp = 0.2_rkind * dtdiff_gp! * this%CFL
        
        if(this%intSharp) then

           if(this%intSharp_gam.lt.-0.5) then !maximize intSharp_gam without restricting time step, based on CFL
              !For intSharp_gam = -1
              !if(this%intSharp_gam.gt.-two) then
              this%mix%intSharp_gam = 0.2_rkind * delta**2 / (dtCFL * this%mix%intSharp_eps + eps)
              !else !ELSE: set intSharp_gam based on maximum velocity as in Tiwari, Freund, Pantano JCP 2013   !For intSharp_gam <= -2
              if(this%intSharp_gam.lt.-1.5) then
                 this%mix%intSharp_gam = zero
                 do i=1,this%mix%ns
                    this%mix%intSharp_gam = max( this%mix%intSharp_gam, P_MAXVAL( four*sqrt( this%u**2 + this%v**2 + this%w**2 ) * this%mix%material(i)%VF*(one-this%mix%material(i)%VF)) )
                 enddo
                 if(this%intSharp_gam.lt.-2.5) then
                    this%mix%intSharp_gam = zero
                    this%mix%intSharp_gam = max( this%mix%intSharp_gam, P_MAXVAL(sqrt( this%u**2 + this%v**2 + this%w**2 ) ) )
                 endif
                 if (this%intSharp_gam.lt.-3.5) then
                 this%mix%intSharp_gam = zero
                 do i=1,this%mix%ns
                     this%mix%intSharp_gam = max( this%mix%intSharp_gam,P_MAXVAL( four*sqrt( this%v**2 ) *this%mix%material(i)%VF*(one-this%mix%material(i)%VF)) )
                 enddo
                 endif

                 if (this%intSharp_gam.lt.-4.5) then
                 this%mix%intSharp_gam = zero
                 do i=1,this%mix%ns
                     this%mix%intSharp_gam = max(this%mix%intSharp_gam,P_MAXVAL( sqrt( this%v**2 )) )
                 enddo
                 endif



              endif
              !print*,this%intSharp_gam,this%mix%intSharp_gam
           endif

           ! if (this%step .LE. this%st_limit) then
           !    !this%mix%intSharp_gam = zero
           !    this%mix%intSharp_gam = this%mix%intSharp_gam * 1.0D-2
           !    if (nrank.eq.0) print*,"limiting intSharp_gam"
           ! endif
           ! ! ! this%mix%intSharp_gam = zero
           ! ! ! if (nrank.eq.0) print*,"limiting intSharp_gam"
              

           dtSharp_diff =  delta**2 / (P_MAXVAL( 6*this%mix%intSharp_gam*this%mix%intSharp_eps) + eps)*this%CFL !based on diffusivity in VF sharpening equation 

           !dtSharp_diff = delta**2/(2.0*this%mix%intSharp_gam*this%mix%intSharp_eps)!from Suhas Jain, Mani, Moin JCP 2020 -- not work

           dtSharp_bound = one/eps
           if(this%intSharp_msk) then
              do i=1,this%mix%ns
                 !dtSharp_bound = min(dtSharp_bound, 0.2_rkind * delta**2 / (P_MAXVAL( this%mix%intSharp_gam*this%mix%intSharp_eps*this%mix%VFboundDiff(:,:,:,i)) + eps))! * this%CFL !based on VF out of bounds diffusivity in VF sharpening equation 
                 !dtSharp_bound = min(dtSharp_bound, 0.2_rkind * delta**2 / (P_MAXVAL( this%mix%intSharp_gam*this%mix%intSharp_eps*this%mix%VFboundDiff(:,:,:,i)) + eps))! * this%CFL !based on VF out of bounds diffusivity in VF sharpening equation 
                 !dtSharp_bound = min(dtSharp_bound, 0.2_rkind * delta**2 / (P_MAXVAL( this%mix%intSharp_dif*this%mix%intSharp_gam*this%intSharp_eps*this%mix%VFboundDiff(:,:,:,i)) + eps))! * this%CFL !based on VF out of bounds diffusivity in VF sharpening equation 
                 dtSharp_bound = min(dtSharp_bound, 0.2_rkind * delta**2 / (P_MAXVAL( this%intSharp_dif*this%mix%intSharp_gam*this%intSharp_eps*this%mix%VFboundDiff(:,:,:,i)) + eps))! * this%CFL !based on VF out of bounds diffusivity in VF sharpening equation 
              enddo
           endif

           dtSharp_Adiff   = delta/max(this%mix%intSharp_gam,eps)! * this%CFL !based on anti-diffusivity in VF sharpening equation
        endif
         a = -1
        ! Use fixed time step if CFL <= 0
        if ( this%CFL .LE. zero ) then
            this%dt = this%dtfixed
            stability = 'fixed'
        else
            stability = 'convective'
            this%dt = dtCFL
            if ( this%dt > dtmu ) then
                 this%dt = dtmu
                 stability = 'shear'
             else if ( this%dt > dtbulk ) then
                 this%dt = dtbulk
                 stability = 'bulk'
             else if ( this%dt > dtYs1 ) then
                 this%dt = dtYs1
                 stability = 'Ys1'
             else if ( this%dt > dtYs2) then
                 this%dt = dtYs2
                 stability = 'Ys2'
             else if ( this%dt > dtVF1 ) then
                 this%dt = dtVF1
                 stability = 'VF1'
             else if ( this%dt > dtVF2) then
                 this%dt = dtVF2
                 stability = 'VF2'
             else if ( this%dt > dtkap ) then
                 this%dt = dtkap
                 stability = 'conductive'
             else if ( this%dt > dtdiff ) then
                 this%dt = dtdiff
                 stability = 'diffusive'
             else if ( this%dt > dtdiff_g ) then
                 this%dt = dtdiff_g
                 stability = 'diffusive g'
             else if ( this%dt > dtdiff_gt ) then
                 this%dt = dtdiff
                 stability = 'diffusive g_t'
             else if ( this%dt > dtplast ) then
                 this%dt = dtplast
                 stability = 'plastic'
             endif
            
            if ( this%dt > dtmu ) then
               this%dt = dtmu
               write(str,'(ES10.3E3)') 1.0D0-dtmu/dtCFL
               stability = 'shear: '//trim(str)//' CFL loss fraction'
            endif
            if ( this%dt > dtbulk ) then
               this%dt = dtbulk
               write(str,'(ES10.3E3)') 1.0D0-dtbulk/dtCFL
               stability = 'bulk: '//trim(str)//' CFL loss fraction'
            endif
            if ( this%dt > dtkap ) then
               this%dt = dtkap
               write(str,'(ES10.3E3)') 1.0D0-dtkap/dtCFL
               stability = 'conductive: '//trim(str)//' CFL loss fraction'
            endif
            if ( this%dt > dtdiff ) then
               this%dt = dtdiff
               write(str,'(ES10.3E3)') 1.0D0-dtdiff/dtCFL
               stability = 'diffusive: '//trim(str)//' CFL loss fraction'
            endif
            if ( this%dt > dtdiff_g ) then
               this%dt = dtdiff_g
               write(str,'(ES10.3E3)') 1.0D0-dtdiff_g/dtCFL
               stability = 'diffusive g: '//trim(str)//' CFL loss fraction'
            endif
            if ( this%dt > dtdiff_gt ) then
               this%dt = dtdiff_gt
               write(str,'(ES10.3E3)') 1.0D0-dtdiff_gt/dtCFL
               stability = 'diffusive g_t: '//trim(str)//' CFL loss fraction'
            endif
            if ( this%dt > dtdiff_gp ) then
               this%dt = dtdiff_gp
               write(str,'(ES10.3E3)') 1.0D0-dtdiff_gp/dtCFL
               stability = 'diffusive g_p: '//trim(str)//' CFL loss fraction'
            endif
	if ((this%use_surfaceTension) .OR. (this%use_CnsrvSurfaceTension)) then	
	    if ( this%dt > dtsigma ) then
               this%dt = dtsigma
               write(str,'(ES10.3E3)') 1.0D0-dtsigma/dtCFL
               stability = 'surfaceTension: '//trim(str)//' CFL loss fraction'
            endif
	endif
            if ( this%dt > dtplast ) then
               this%dt = dtplast
               write(str,'(ES10.3E3)') 1.0D0-dtplast/dtCFL
               stability = 'plastic: '//trim(str)//' CFL loss fraction'
            end if
            if (this%intSharp) then
               if ( this%dt > dtSharp_diff ) then
                  this%dt = dtSharp_diff
                  write(str,'(ES10.3E3)') 1.0D0-dtSharp_diff/dtCFL
                  stability = 'sharp diff: '//trim(str)//' CFL loss fraction'
               end if
               if ( this%dt > dtSharp_Adiff ) then
                  this%dt = dtSharp_Adiff
                  write(str,'(ES10.3E3)') 1.0D0-dtSharp_Adiff/dtCFL
                  stability = 'sharp a-diff: '//trim(str)//' CFL loss fraction'
               end if
               
               if(this%intSharp_msk) then

                  if ( this%dt > dtSharp_bound ) then
                      this%dt = dtSharp_bound
                      !write(str2,'(F25.18)') dtCFL
                      !write(str,'(F6.2)') 1.0D2*dtSharp_bound/dtCFL
                      !stability = 'Sharp VF bounds: '//trim(str)//'%'//' '//trim(str2)
                      write(str,'(ES10.3E3)') 1.0D0-dtSharp_bound/dtCFL
                      stability = 'sharp VF bounds: '//trim(str)//' CFL loss fraction'
                  end if
               end if
            end if

            if (this%step .LE. this%st_limit) then
               this%dt = min(this%dt / st_fac, this%dtfixed)
               stability = 'startup'
            endif
         endif

    end subroutine

    subroutine get_primitive(this)
        use reductions, only: P_MAXVAL, P_MINVAL
        use decomp_2d,  only: nrank
        class(sgrid), target, intent(inout) :: this
        real(rkind), dimension(:,:,:), pointer :: onebyrho
       ! real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: onebyrho
        real(rkind), dimension(:,:,:), pointer :: rhou,rhov,rhow,TE
        real(rkind) :: rhomin

        onebyrho => this%ybuf(:,:,:,1)

        call this%mix%get_rho(this%rho)

        ! rhomin = P_MINVAL(this%rho)
        ! if (nrank.eq.0) print*,rhomin

        rhou => this%Wcnsrv(:,:,:,mom_index  )
        rhov => this%Wcnsrv(:,:,:,mom_index+1)
        rhow => this%Wcnsrv(:,:,:,mom_index+2)
        TE   => this%Wcnsrv(:,:,:, TE_index  )

        onebyrho = one/this%rho
        this%u = rhou * onebyrho
        this%v = rhov * onebyrho
        this%w = rhow * onebyrho

        if(this%use_CnsrvSurfaceTension) then
     
           this%e = ((TE-this%mix%surfaceTension_pe)*onebyrho) - half*( this%u*this%u + this%v*this%v + this%w*this%w )

        else
           
           this%e = (TE - 0.5*(rhou*rhou +rhov*rhov + rhow*rhow)/this%rho)/this%rho !(TE*onebyrho) - half*(this%u*this%u + this%v*this%v +this%w*this%w )

        endif
    
        call this%mix%get_primitive(this%rho, this%u, this%v, this%w, this%e, this%devstress, this%p, this%sos)                  ! Get primitive variables for individual species

    end subroutine


    subroutine get_primitive_g(this)
      class(sgrid), target, intent(inout) :: this

      call this%mix%get_primitive_g(this%rho)                  ! Get primitive kinematic variables
      
    end subroutine get_primitive_g
     

    pure subroutine get_conserved(this)
        class(sgrid), intent(inout) :: this

        ! Assume rho is already available
        this%Wcnsrv(:,:,:,mom_index  ) = this%rho * this%u
        this%Wcnsrv(:,:,:,mom_index+1) = this%rho * this%v
        this%Wcnsrv(:,:,:,mom_index+2) = this%rho * this%w

        if(this%use_CnsrvSurfaceTension) then
            this%Wcnsrv(:,:,:, TE_index  ) = this%rho * ( this%e + half*(this%u*this%u + this%v*this%v + this%w*this%w )) + this%mix%surfaceTension_pe
        else
            this%Wcnsrv(:,:,:, TE_index  ) = this%rho * ( this%e +half*(this%u*this%u + this%v*this%v + this%w*this%w ))
        endif
        ! add 2M (mass fraction and hydrodynamic energy) variables here
        call this%mix%get_conserved(this%rho,this%u,this%v,this%w)

    end subroutine

    subroutine get_conserved_g(this)
      class(sgrid), target, intent(inout) :: this

      call this%mix%get_conserved_g(this%rho)                  ! Get conserved kinematic variables
      
    end subroutine get_conserved_g

    subroutine post_bc(this)
        class(sgrid), intent(inout) :: this

        if(this%useOneG) then
            call this%mix%get_mixture_properties()
        endif
        call this%mix%get_eelastic_devstress(this%devstress)   ! Get species elastic energies, and mixture and species devstress
        call this%mix%get_ehydro_from_p(this%rho)              ! Get species hydrodynamic energy, temperature; and mixture pressure, temperature
        call this%mix%get_pmix(this%p)                         ! Get mixture pressure
        call this%mix%get_Tmix(this%T)                         ! Get mixture temperature
        call this%mix%getSOS(this%rho,this%p,this%sos)
!print *, 'SOS: ', this%sos(179,1,1)
        ! assuming pressures have relaxed and sum( (Ys*(ehydro + eelastic) ) over all
        ! materials equals e
        call this%mix%get_emix(this%e)
    end subroutine

    subroutine post_bc_2(this)
        class(sgrid), intent(inout) :: this

        if(this%useOneG) then
            call this%mix%get_mixture_properties()
        endif
        call this%mix%get_eelastic_devstress(this%devstress)   ! Get specieselastic energies, and mixture and species devstress
        ! Get specieshydrodynamic energy, temperature; and mixture pressure, temperature
        call this%mix%get_ehydro_from_p(this%rho) 
        call this%mix%get_pmix(this%p)                         ! Get mixturepressure
        call this%mix%get_Tmix(this%T)                         ! Get mixturetemperature
        call this%mix%getSOS(this%rho,this%p,this%sos)
!print *, 'SOS: ', this%sos(179,1,1)
        ! assuming pressures have relaxed and sum( (Ys*(ehydro + eelastic) )
        ! over all
        ! materials equals e
        call this%mix%get_emix(this%e)

    end subroutine

    subroutine CheckTau(this,tauxx,tauyy,tauzz,tauxy,tauyx,tauyz,tauzy,tauxz,tauzx)
        use decomp_2d, only: transpose_y_to_x, transpose_x_to_y,transpose_y_to_z, transpose_z_to_y
        use exits,      only: message,nancheck,GracefulExit
        class(sgrid), target, intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in) :: tauxx,tauyy,tauzz,tauxy,tauyx,tauyz,tauzy,tauxz,tauzx
        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: x_half, y_half,z_half,tau11, tau22, tau33, tau12, tau21, tau13, tau31, tau23, tau32

        print *, "Check Tau"
        x_half = this%x + 0.5*this%dx 
        y_half = this%y + 0.5*this%dy
        z_half = this%z + 0.5*this%dz
          
        tau11 = -2_rkind / 3._rkind*this%mu*(-4.0*sin(4.0*this%z)*cos(x_half)) + 4._rkind / 3._rkind*this%mu*4.0*sin(2.0*this%z)*cos(4.0*x_half)
        tau33 = 4._rkind / 3._rkind*this%mu*(-4.0*sin(4.0*z_half)*cos(this%x)) - 2._rkind /3._rkind*this%mu*(4.0*sin(2.0*z_half)*cos(4.0*this%x)) 
        tau22 = 0; tau12 = 0; tau21 = 0; tau23 = 0; tau32 = 0;
        tau13 = this%mu*(2*cos(2*z_half)*sin(4*this%x) - cos(4*z_half)*sin(this%x))
        tau31 = this%mu*(2*cos(2*this%z)*sin(4*x_half) - cos(4*this%z)*sin(x_half)) 
        
        this%tauxx = tauxx; this%tauxxe = tau11; 
        this%tauxy = tauxy; this%tauxye = tau12;
        this%tauyx = tauyx; this%tauyxe = tau21;
        this%tauyy = tauyy; this%tauyye = tau22; 
        this%tauyz = tauyz; this%tauyze = tau23;
        this%tauzy = tauzy; this%tauzye = tau32;
        this%tauxz = tauxz; this%tauxze = tau13;
        this%tauzx = tauzx; this%tauzxe = tau31;
    end subroutine

    subroutine getRHS(this, rhs, divu, viscwork)
        use decomp_2d, only: transpose_y_to_x, transpose_x_to_y,transpose_y_to_z, transpose_z_to_y
        use operators, only: divergence,gradient,divergenceFV, interpolateFV, interpolateFV_x, interpolateFV_y, interpolateFV_z, gradFV_N2Fx, gradFV_N2Fy, gradFV_N2Fz
        use exits,      only: message,nancheck,GracefulExit
        class(sgrid), target, intent(inout) :: this
        real(rkind), dimension(this%nxp, this%nyp, this%nzp,ncnsrv), intent(out) :: rhs
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),     intent(out) :: divu
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),target, intent(out) :: viscwork
        real(rkind), dimension(this%nxp, this%nyp, this%nzp,9), target :: duidxj, duidxj_s
        real(rkind), dimension(this%nxp, this%nyp, this%nzp,12), target :: duidxj_int        
        real(rkind), dimension(:,:,:), pointer :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
        real(rkind), dimension(:,:,:), pointer :: dudx_s,dudy_s,dudz_s,dvdx_s,dvdy_s,dvdz_s,dwdx_s,dwdy_s,dwdz_s
        real(rkind), dimension(:,:,:), pointer :: dvdy_x,dwdz_x, dvdx_y, dwdx_z
        real(rkind), dimension(:,:,:), pointer :: dudx_y, dwdz_y, dudy_x, dwdy_z
        real(rkind), dimension(:,:,:), pointer :: dudx_z, dvdy_z, dudz_x, dvdz_y
        real(rkind), dimension(:,:,:), pointer :: tauxx,tauxy,tauxz,tauyy,tauyz,tauzz, tauzy, tauzx, tauyx
        real(rkind), dimension(:,:,:), pointer :: qx,qy,qz
        real(rkind), dimension(:,:,:), pointer :: ehmix
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,3) :: u_int, v_int,w_int
        integer :: imat,i
        !logical :: useNewSPF = .FALSE.
        real(rkind), dimension(this%nxp, this%nyp, this%nzp) :: ke,tmp,dJ,drhodx,drhody,drhodz, uJ, vJ, wJ, keJ, eJ, Fbody
        real(rkind), dimension(this%nxp, this%nyp, this%nzp) :: drhoedx,drhoedy, drhoedz
        real(rkind), dimension(this%nxp, this%nyp, this%nzp, 3) :: J,Frho,Fenergy, ke_int, Fp, yMetric_F2N_int
        real(rkind) :: g = -0.1

        !this%u = sin(2*this%y)*sin(4*this%x)
        !this%v = 0
        !this%w = cos(4*this%y)*cos(this%x)
        !sin(2*this%y)*sin(4*this%x); this%v = 0; this%w =  0; !cos(4*this%y)*cos(this%x); this%w = 0;

        if(this%use_Stagg) then

           dudx => duidxj(:,:,:,1); dudy => duidxj(:,:,:,2); dudz => duidxj(:,:,:,3);
           dvdx => duidxj(:,:,:,4); dvdy => duidxj(:,:,:,5); dvdz => duidxj(:,:,:,6);
           dwdx => duidxj(:,:,:,7); dwdy => duidxj(:,:,:,8); dwdz => duidxj(:,:,:,9);

           dvdy_x => duidxj_int(:,:,:,1); dudx_y => duidxj_int(:,:,:,2); dudx_z => duidxj_int(:,:,:,3);
           dwdz_x => duidxj_int(:,:,:,4); dwdz_y => duidxj_int(:,:,:,5); dvdy_z => duidxj_int(:,:,:,6);
           dvdx_y => duidxj_int(:,:,:,7); dudy_x => duidxj_int(:,:,:,8); dudz_x => duidxj_int(:,:,:,9);
           dwdx_z => duidxj_int(:,:,:,10); dwdy_z => duidxj_int(:,:,:,11); dvdz_y => duidxj_int(:,:,:,12);

           
           call gradient(this%decomp,this%derCD06_nostretch,this%u, dudx, dudy, dudz,  -this%x_bc,  this%y_bc,this%z_bc)
           call gradient(this%decomp,this%derCD06_nostretch,this%v, dvdx, dvdy, dvdz,  this%x_bc, -this%y_bc,this%z_bc)
           call gradient(this%decomp,this%derCD06_nostretch,this%w, dwdx, dwdy, dwdz,  this%x_bc,  this%y_bc,-this%z_bc)

           call interpolateFV_x(this%decomp,this%interpMid,dvdy,dvdy_x,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
           call interpolateFV_x(this%decomp,this%interpMid,dwdz,dwdz_x,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
           call interpolateFV_y(this%decomp,this%interpMid,dvdx,dvdx_y,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
           call interpolateFV_z(this%decomp,this%interpMid,dwdx,dwdx_z,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)

           call interpolateFV_y(this%decomp,this%interpMid,dudx,dudx_y,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
           call interpolateFV_y(this%decomp,this%interpMid,dwdz,dwdz_y,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
           call interpolateFV_x(this%decomp,this%interpMid,dudy,dudy_x,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
           call interpolateFV_z(this%decomp,this%interpMid,dwdy,dwdy_z,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)

           call interpolateFV_z(this%decomp,this%interpMid,dudx,dudx_z,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
           call interpolateFV_z(this%decomp,this%interpMid,dvdy,dvdy_z,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc) 
           call interpolateFV_x(this%decomp,this%interpMid,dudz,dudz_x,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
           call interpolateFV_y(this%decomp,this%interpMid,dvdz,dvdz_y,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)

           if(this%Stretch1Dy) then

              call interpolateFV(this%decomp,this%interpMid,this%yMetric_F2N,yMetric_F2N_int,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
              dudy_x = yMetric_F2N_int(:,:,:,1)*dudy_x
              dvdy_x = yMetric_F2N_int(:,:,:,1)*dvdy_x
              dwdy_z = yMetric_F2N_int(:,:,:,3)*dwdy_z
              dvdy_z = yMetric_F2N_int(:,:,:,3)*dvdy_z 
           endif
         
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! GET STAGGERED DERIVATIVES !!!!!!!!!!!!!!!!!!!!
           dudx_s => duidxj_s(:,:,:,1); dudy_s => duidxj_s(:,:,:,2); dudz_s => duidxj_s(:,:,:,3);
           dvdx_s => duidxj_s(:,:,:,4); dvdy_s => duidxj_s(:,:,:,5); dvdz_s => duidxj_s(:,:,:,6);
           dwdx_s => duidxj_s(:,:,:,7); dwdy_s => duidxj_s(:,:,:,8); dwdz_s => duidxj_s(:,:,:,9);

           call gradFV_N2Fx(this%decomp,this%derStagg,this%u,dudx_s,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
           call gradFV_N2Fx(this%decomp,this%derStagg,this%v,dvdx_s,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
           call gradFV_N2Fx(this%decomp,this%derStagg,this%w,dwdx_s,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)

           call gradFV_N2Fy(this%decomp,this%derStagg,this%u,dudy_s,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
           call gradFV_N2Fy(this%decomp,this%derStagg,this%v,dvdy_s,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
           call gradFV_N2Fy(this%decomp,this%derStagg,this%w,dwdy_s,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
          
           call gradFV_N2Fz(this%decomp,this%derStagg,this%u,dudz_s,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
           call gradFV_N2Fz(this%decomp,this%derStagg,this%v,dvdz_s,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
           call gradFV_N2Fz(this%decomp,this%derStagg,this%w,dwdz_s,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc) 


        else

           dudx => duidxj(:,:,:,1); dudy => duidxj(:,:,:,2); dudz => duidxj(:,:,:,3);
           dvdx => duidxj(:,:,:,4); dvdy => duidxj(:,:,:,5); dvdz => duidxj(:,:,:,6);
           dwdx => duidxj(:,:,:,7); dwdy => duidxj(:,:,:,8); dwdz => duidxj(:,:,:,9);

           call this%gradient(this%u, dudx, dudy, dudz, -this%x_bc,  this%y_bc,  this%z_bc)
           call this%gradient(this%v, dvdx, dvdy, dvdz,  this%x_bc, -this%y_bc,  this%z_bc)
           call this%gradient(this%w, dwdx, dwdy, dwdz,  this%x_bc,  this%y_bc, -this%z_bc)
        endif

        divu = dudx + dvdy + dwdz

        call this%getPhysicalProperties()
        !call this%LAD%get_viscosities(this%rho,duidxj,this%mu,this%bulk,this%x_bc,this%y_bc,this%z_bc)
        call this%LAD%get_viscosities(this%rho,this%p,this%sos,duidxj,this%mu,this%bulk,this%x_bc,this%y_bc,this%z_bc,this%dt,this%intSharp_pfloor,this%yMetric,this%dy_stretch)
        call this%LAD%get_conductivity(this%rho,this%p,this%e,this%T,this%sos,this%kap,this%x_bc,this%y_bc,this%z_bc,this%intSharp_tfloor)
        !call this%LAD%get_P_conductivity(this%rho,this%p,this%e,this%T,this%sos,this%kap,this%x_bc,this%y_bc,this%z_bc,this%intSharp_tfloor)
        if (this%PTeqb) then
            ! subtract elastic energies to determine mixture hydrostatic energy. conductivity 
            ! is assumed a function of only hydrostatic energy
            ehmix => viscwork ! use some storage space
            ehmix = this%e
            do imat = 1, this%mix%ns
                ehmix = ehmix - this%mix%material(imat)%Ys * this%mix%material(imat)%eel
            enddo
            call this%LAD%get_conductivity(this%rho,this%p,ehmix,this%T,this%sos,this%kap,this%x_bc,this%y_bc,this%z_bc,this%intSharp_tfloor)
        end if

        if( .NOT. this%use_Stagg) then

           ! Get tau tensor tensor. Put in off-diagonal components of duidxj (also get the viscous work term for energy equation)
           call this%get_tau( duidxj, viscwork )
           ! Now, associate the pointers to understand what's going on better
           tauxx => duidxj(:,:,:,tauxxidx); tauxy => duidxj(:,:,:,tauxyidx); tauxz => duidxj(:,:,:,tauxzidx);
                                         tauyy => duidxj(:,:,:,tauyyidx); tauyz => duidxj(:,:,:,tauyzidx);
                                                                          tauzz => duidxj(:,:,:,tauzzidx);
           !print '(a,9(e21.14,1x))', 'dudx ', duidxj(89,1,1,1:9)
           !print '(a,9(e21.14,1x))', 'tauxx', tauxx(89,1,1)
           ! Add the deviatoric stress to the tau for use in fluxes 
           tauxx = tauxx + this%sxx; tauxy = tauxy + this%sxy; tauxz = tauxz + this%sxz
                                  tauyy = tauyy + this%syy; tauyz = tauyz + this%syz
                                                            tauzz = tauzz + this%szz
           !print '(a,9(e21.14,1x))', 'tauxx', tauxx(89,1,1)
      
           ! store artificial stress tensor in devstress. this should not break anything since devstress will be
           ! overwritten in get_primitive and post_bc. used in update_eh -- NSG
           this%sxx = tauxx - this%sxx; this%sxy = tauxy - this%sxy; this%sxz = tauxz - this%sxz
                                     this%syy = tauyy - this%syy; this%syz = tauyz - this%syz
                                                                  this%szz = tauzz - this%szz
      
           ! Get heat conduction vector (q). Stored in remaining 3 components of duidxj 
           qx => duidxj(:,:,:,qxidx); qy => duidxj(:,:,:,qyidx); qz => duidxj(:,:,:,qzidx);

        else
           call this%get_tauStagg( duidxj, duidxj_int, duidxj_s )    
           tauxx => duidxj(:,:,:,tauxxidx); tauyx => duidxj(:,:,:,tauxyidx); tauzx => duidxj(:,:,:,tauxzidx);
           tauxy => duidxj_int(:,:,:, 1);   tauyy => duidxj(:,:,:,tauyyidx); tauzy => duidxj(:,:,:,tauyzidx);
           tauxz => duidxj_int(:,:,:,2);    tauyz => duidxj_int(:,:,:,3);    tauzz => duidxj(:,:,:,tauzzidx);
           ! Add the deviatoric stress to the tau for use in fluxes 
           !tauxx = tauxx + this%sxx; tauxy = tauxy + this%sxy; tauxz = tauxz + this%sxz
           !                       tauyy = tauyy + this%syy; tauyz = tauyz + this%syz
           !                                                 tauzz = tauzz + this%szz
          
           ! store artificial stress tensor in devstress. this should not break
           ! anything since devstress will be
           ! overwritten in get_primitive and post_bc. used in update_eh -- NSG
           !this%sxx = tauxx - this%sxx; this%sxy = tauxy - this%sxy; this%sxz = tauxz - this%sxz
           !                            this%syy = tauyy - this%syy; this%syz = tauyz - this%syz
           !                                                          this%szz = tauzz - this%szz
           this%tauxx = tauxx; this%tauyy = tauyy; this%tauzz = tauzz;
           this%tauxy = tauxy; this%tauyx = tauyx; this%tauxz = tauxz;
           this%tauzx = tauzx; this%tauyz = tauyz; this%tauzy = tauzy;
           ! Get heat conduction vector (q). Stored in remaining 3 components of
           ! duidxj 
           qx => duidxj(:,:,:,qxidx); qy => duidxj(:,:,:,qyidx); qz => duidxj(:,:,:,qzidx);

           !call this%CheckTau(tauxx,tauyy,tauzz,tauxy,tauyx,tauyz,tauzy,tauxz,tauzx)

        endif 

      !  call this%mix%get_qmix(qx, qy, qz)                     ! Get only species diffusion fluxes if PTeqb, else, everything
        if (this%PTeqb) then
          call this%get_q(qx, qy, qz)            ! add artificial thermal conduction fluxes
        end if

        call this%get_qLAD(qx,qy,qz,this%qDiv)
        rhs = zero

        if(this%use_Stagg) then

              call this%getRHS_xStagg(              rhs,&
                                      tauxx,tauxy,tauxz,&
                                                     qx )

              call this%getRHS_yStagg(              rhs,&
                                      tauxy,tauyy,tauyz,&
                                                     qy )
              call this%getRHS_zStagg(              rhs,&
                                      tauxz,tauyz,tauzz,&
                                                     qz )




           else
              call this%getRHS_x(              rhs,&
                                 tauxx,tauxy,tauxz,&
                                                qx )
           !print '(a,4(e21.14,1x))', 'rhsx: ', rhs(179,1,1,1:4)
              call this%getRHS_y(              rhs,&
                                 tauxy,tauyy,tauyz,&
                                                qy )
           !print '(a,4(e21.14,1x))', 'rhsy: ', rhs(179,1,1,1:4)

              call this%getRHS_z(              rhs,&
                                 tauxz,tauyz,tauzz,&
                                                 qz )
           !print '(a,4(e21.14,1x))', 'rhsz: ', rhs(179,1,1,1:4)
        endif

        if(this%intSharp .AND. this%intSharp_cpl) then
           !calculate kinetic energy for intSharp terms
            ke = half*( this%u**2 + this%v**2 + this%w**2 ) 
       
             !call divergence(this%decomp,this%derD02,this%mix%intSharp_f(:,:,:,1)*this%u,this%mix%intSharp_f(:,:,:,2)*this%u,this%mix%intSharp_f(:,:,:,3)*this%u,tmp,this%x_bc,-this%y_bc,-this%z_bc)
             ! rhs(:,:,:,mom_index  ) = rhs(:,:,:,mom_index  ) + tmp

             !call divergence(this%decomp,this%derD02,this%mix%intSharp_f(:,:,:,1)*this%v,this%mix%intSharp_f(:,:,:,2)*this%v,this%mix%intSharp_f(:,:,:,3)*this%v,tmp,-this%x_bc,this%y_bc,-this%z_bc)
             ! rhs(:,:,:,mom_index+1) = rhs(:,:,:,mom_index+1) + tmp

             ! call divergence(this%decomp,this%derD02,this%mix%intSharp_f(:,:,:,1)*this%w,this%mix%intSharp_f(:,:,:,2)*this%w,this%mix%intSharp_f(:,:,:,3)*this%w,tmp,-this%x_bc,-this%y_bc,this%z_bc)
             ! rhs(:,:,:,mom_index+2) = rhs(:,:,:,mom_index+2) + tmp

             ! call divergence(this%decomp,this%derD02,this%mix%intSharp_f(:,:,:,1)*ke +this%mix%intSharp_h(:,:,:,1),this%mix%intSharp_f(:,:,:,2)*ke +this%mix%intSharp_h(:,:,:,2),this%mix%intSharp_f(:,:,:,3)*ke +this%mix%intSharp_h(:,:,:,3),tmp,-this%x_bc,-this%y_bc,-this%z_bc)
              !rhs(:,:,:,TE_index   ) = rhs(:,:,:,TE_index   ) + tmp


              !high order terms
              !call divergence(this%decomp,this%der,this%mix%intSharp_fDiff(:,:,:,1)*this%u,this%mix%intSharp_fDiff(:,:,:,2)*this%u,this%mix%intSharp_fDiff(:,:,:,3)*this%u,tmp,this%x_bc,-this%y_bc,-this%z_bc)
              !rhs(:,:,:,mom_index  ) = rhs(:,:,:,mom_index  ) + tmp

              !call divergence(this%decomp,this%der,this%mix%intSharp_fDiff(:,:,:,1)*this%v,this%mix%intSharp_fDiff(:,:,:,2)*this%v,this%mix%intSharp_fDiff(:,:,:,3)*this%v,tmp,-this%x_bc,this%y_bc,-this%z_bc)
              !rhs(:,:,:,mom_index+1) = rhs(:,:,:,mom_index+1) + tmp

              !call divergence(this%decomp,this%der,this%mix%intSharp_fDiff(:,:,:,1)*this%w,this%mix%intSharp_fDiff(:,:,:,2)*this%w,this%mix%intSharp_fDiff(:,:,:,3)*this%w,tmp,-this%x_bc,-this%y_bc,this%z_bc)
              !rhs(:,:,:,mom_index+2) = rhs(:,:,:,mom_index+2) + tmp

              !call divergence(this%decomp,this%der,this%mix%intSharp_fDiff(:,:,:,1)*ke +this%mix%intSharp_hDiff(:,:,:,1),this%mix%intSharp_fDiff(:,:,:,2)*ke +this%mix%intSharp_hDiff(:,:,:,2),this%mix%intSharp_fDiff(:,:,:,3)*ke +this%mix%intSharp_hDiff(:,:,:,3),tmp,-this%x_bc,-this%y_bc,-this%z_bc)
              !rhs(:,:,:,TE_index   ) = rhs(:,:,:,TE_index   ) + tmp
 



          !FV sharpening
           rhs(:,:,:,mom_index  ) = rhs(:,:,:,mom_index  ) + this%mix%intSharp_fFV(:,:,:,1) + this%mix%intSharp_fDiffFV(:,:,:,1)
           rhs(:,:,:,mom_index+1) = rhs(:,:,:,mom_index+1) + this%mix%intSharp_fFV(:,:,:,2) + this%mix%intSharp_fDiffFV(:,:,:,2)
           rhs(:,:,:,mom_index+2) = rhs(:,:,:,mom_index+2) + this%mix%intSharp_fFV(:,:,:,3) + this%mix%intSharp_fDiffFV(:,:,:,3)
           rhs(:,:,:,TE_index   ) = rhs(:,:,:,TE_index   ) + this%mix%intSharp_hFV + this%mix%intSharp_kFV + this%mix%intSharp_hDiffFV + this%mix%intSharp_kDiffFV

        endif

        if (this%use_surfaceTension .AND. (.NOT. this%use_CnsrvSurfaceTension)) then
            rhs(:,:,:,mom_index  ) = rhs(:,:,:,mom_index  ) + this%mix%surfaceTension_f(:,:,:,1)
            rhs(:,:,:,mom_index+1) = rhs(:,:,:,mom_index+1) + this%mix%surfaceTension_f(:,:,:,2)
            rhs(:,:,:,mom_index+2) = rhs(:,:,:,mom_index+2) + this%mix%surfaceTension_f(:,:,:,3)
            rhs(:,:,:,TE_index   ) = rhs(:,:,:,TE_index   ) + this%mix%surfaceTension_e
        endif

        this%uJ = 0.0; this%vJ = 0.0; this%wJ = 0.0; this%keJ = 0.0; this%eJ = 0; 
        if( .NOT. this%twoPhaseLAD) then
             
          J = 0
          ke = half*( this%u**2 + this%v**2 + this%w**2 )

          do i = 1,this%mix%ns
              J = J + this%mix%material(i)%Ji
          enddo

          call divergence(this%decomp, this%der,J(:,:,:,1), J(:,:,:,2), J(:,:,:,3),dJ,this%x_bc,this%y_bc,this%z_bc )
          uJ = this%u*dJ
          vJ = this%v*dJ
          wJ = this%w*dJ
          keJ = ke*dJ

        else
         
          ke = half*( this%u**2 + this%v**2 + this%w**2 )
          call this%mix%getLAD_5eqn(this%rho,this%p,this%e,Frho,Fenergy,Fp,this%x_bc,this%y_bc,this%z_bc,this%dx,this%dy,this%dz,this%periodicx,this%periodicy,this%periodicz)

          
          if( this%LADInt .OR. this%LADN2F) then
             call interpolateFV(this%decomp,this%interpMid,this%u,u_int,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
             call interpolateFV(this%decomp,this%interpMid,this%v,v_int,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
             call interpolateFV(this%decomp,this%interpMid,this%w,w_int,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)

             ke_int = half*(u_int**2 + v_int**2 + w_int**2)
             call divergenceFV(this%decomp,this%derStagg,u_int(:,:,:,1)*Frho(:,:,:,1),u_int(:,:,:,2)*Frho(:,:,:,2),u_int(:,:,:,2)*Frho(:,:,:,3),this%uJ,this%periodicx, this%periodicy, this%periodicz,this%x_bc,this%y_bc,this%z_bc)
             call divergenceFV(this%decomp,this%derStagg,v_int(:,:,:,1)*Frho(:,:,:,1),v_int(:,:,:,2)*Frho(:,:,:,2),v_int(:,:,:,3)*Frho(:,:,:,3),this%vJ,this%periodicx, this%periodicy, this%periodicz,this%x_bc,this%y_bc,this%z_bc)
             call divergenceFV(this%decomp,this%derStagg,w_int(:,:,:,1)*Frho(:,:,:,1),w_int(:,:,:,2)*Frho(:,:,:,2),w_int(:,:,:,3)*Frho(:,:,:,3),this%wJ,this%periodicx, this%periodicy, this%periodicz,this%x_bc,this%y_bc,this%z_bc)
             call divergenceFV(this%decomp,this%derStagg,ke_int(:,:,:,1)*Frho(:,:,:,1),ke_int(:,:,:,2)*Frho(:,:,:,2),ke_int(:,:,:,3)*Frho(:,:,:,3),this%keJ,this%periodicx, this%periodicy, this%periodicz,this%x_bc,this%y_bc,this%z_bc)
             call divergenceFV(this%decomp,this%derStagg,Fenergy(:,:,:,1),Fenergy(:,:,:,2),Fenergy(:,:,:,3),this%eJ,this%periodicx, this%periodicy, this%periodicz,this%x_bc,this%y_bc,this%z_bc)
             call divergenceFV(this%decomp,this%derStagg,Fp(:,:,:,1),Fp(:,:,:,2),Fp(:,:,:,3),this%pJ,this%periodicx,this%periodicy, this%periodicz,this%x_bc,this%y_bc,this%z_bc)

           else

            call divergence(this%decomp,this%der,this%u*Frho(:,:,:,1),this%u*Frho(:,:,:,2),this%u*Frho(:,:,:,3),this%uJ,this%x_bc,this%y_bc,this%z_bc)
            call divergence(this%decomp,this%der,this%v*Frho(:,:,:,1),this%v*Frho(:,:,:,2),this%v*Frho(:,:,:,3),this%vJ,this%x_bc,this%y_bc,this%z_bc)
            call divergence(this%decomp,this%der,this%w*Frho(:,:,:,1),this%w*Frho(:,:,:,2),this%w*Frho(:,:,:,3),this%wJ,this%x_bc,this%y_bc,this%z_bc)
            call divergence(this%decomp,this%der,ke*Frho(:,:,:,1),ke*Frho(:,:,:,2),ke*Frho(:,:,:,3),this%keJ,this%x_bc,this%y_bc,this%z_bc)
            call divergence(this%decomp,this%der,Fenergy(:,:,:,1),Fenergy(:,:,:,2),Fenergy(:,:,:,3),this%eJ,this%x_bc,this%y_bc,this%z_bc)

           endif 
          rhs(:,:,:, mom_index   ) = rhs(:,:,:,mom_index   ) + this%uJ
          rhs(:,:,:, mom_index+1 ) = rhs(:,:,:,mom_index+1 ) + this%vJ
          rhs(:,:,:, mom_index+2 ) = rhs(:,:,:,mom_index+2 ) + this%wJ
          rhs(:,:,:, TE_index )    = rhs(:,:,:,TE_index    ) + this%keJ + this%eJ

        endif
        !rhs(:,:,:, mom_index+1 ) = rhs(:,:,:,mom_index+1 ) + this%rho*g
        !rhs(:,:,:, TE_index )    = rhs(:,:,:,TE_index    ) + this%v*this%rho*g
        ! Call problem source hook
        call hook_mixture_source(this%decomp, this%mesh, this%fields, this%mix, this%tsim, rhs)
 
    end subroutine



subroutine getRHS_NC(this, rhs, divu, viscwork)
        use operators, only: divergence,gradient
        use exits, only: GracefulExit
        class(sgrid), target, intent(inout) :: this
        real(rkind), dimension(this%nxp, this%nyp, this%nzp,ncnsrv), intent(out):: rhs
        real(rkind), dimension(this%nxp, this%nyp, this%nzp,9), target :: duidxj 
        real(rkind), dimension(this%nxp, this%nyp, this%nzp,3*this%mix%ns),target :: gradYs
        real(rkind), dimension(this%nxp, this%nyp, this%nzp,9), target :: NCbuff
! extra buffers for nonconservative
        real(rkind), dimension(this%nxp, this%nyp, this%nzp) :: flux, tmp, ke
        real(rkind), dimension(:,:,:), pointer :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
        real(rkind), dimension(:,:,:), pointer :: d2udx2,d2udy2,d2udz2,d2vdx2,d2vdy2,d2vdz2,d2wdx2,d2wdy2,d2wdz2
        real(rkind), dimension(:,:,:,:), pointer :: dYsdx, dYsdy, dYsdz
        real(rkind), dimension(:,:,:,:), pointer :: d2Ysdx2, d2Ysdy2, d2Ysdz2
        real(rkind), dimension(:,:,:), pointer :: tauxx,tauxy,tauxz,tauyy,tauyz,tauzz
        real(rkind), dimension(:,:,:), pointer :: dTdx,dTdy,dTdz,d2Tdx2,d2Tdy2,d2Tdz2
        real(rkind), dimension(:,:,:), pointer :: drhou_dx, drhov_dy, drhow_dz
        real(rkind), dimension(:,:,:), pointer :: dmudx,dmudy,dmudz
        real(rkind), dimension(:,:,:), pointer :: totale,bambda,lambda
        real(rkind), dimension(:,:,:), pointer :: dbulkdx,dbulkdy,dbulkdz
        real(rkind), dimension(:,:,:), pointer :: dkapdx,dkapdy,dkapdz
        real(rkind), dimension(:,:,:,:), pointer :: Jx,Jy,Jz
        real(rkind), dimension(:,:,:), pointer :: dsumJxdx, dsumJydy, dsumJzdz
        real(rkind), dimension(:,:,:), pointer :: dJxdx, dJydy, dJzdz
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),target, intent(out) :: viscwork
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),     intent(out) :: divu
        real(rkind), dimension(:,:,:), pointer :: ehmix
        real(rkind), dimension(:,:,:), pointer :: qx,qy,qz
        integer :: i,k, imat
        
        type(derivatives), pointer :: der
        type(decomp_info), pointer :: decomp
        real(rkind), dimension(:,:,:), pointer :: xtmp1,xtmp2,ztmp1,ztmp2

        der => this%der
        decomp => this%decomp
        xtmp1 => this%xbuf(:,:,:,1)
        xtmp2 => this%xbuf(:,:,:,2)
        ztmp1 => this%zbuf(:,:,:,1)
        ztmp2 => this%zbuf(:,:,:,2)

        dudx => duidxj(:,:,:,1); dudy => duidxj(:,:,:,2); dudz => duidxj(:,:,:,3);
        dvdx => duidxj(:,:,:,4); dvdy => duidxj(:,:,:,5); dvdz => duidxj(:,:,:,6);
        dwdx => duidxj(:,:,:,7); dwdy => duidxj(:,:,:,8); dwdz => duidxj(:,:,:,9);

        call this%gradient(this%u,dudx,dudy,dudz, [0,0], this%y_bc, this%z_bc)
        call this%gradient(this%v,dvdx,dvdy,dvdz, this%x_bc, [0,0], this%z_bc)
        call this%gradient(this%w,dwdx,dwdy,dwdz, this%x_bc, this%y_bc, [0,0])

        divu = dudx + dvdy + dwdz




        if (this%mix%ns .GT. 1) then
          dYsdx => gradYs(:,:,:,1:this%mix%ns); dYsdy => gradYs(:,:,:,this%mix%ns+1:2*this%mix%ns);
          dYsdz => gradYs(:,:,:,2*this%mix%ns+1:3*this%mix%ns);
            do i = 1,this%mix%ns
                call this%gradient(this%mix%material(i)%Ys,dYsdx(:,:,:,i),dYsdy(:,:,:,i),dYsdz(:,:,:,i),this%x_bc, this%y_bc, this%z_bc)
         end do
       endif

        call this%getPhysicalProperties()
        !call
        !this%LAD%get_viscosities(this%rho,duidxj,this%mu,this%bulk,this%x_bc,this%y_bc,this%z_bc)
        call this%LAD%get_viscosities(this%rho,this%p,this%sos,duidxj,this%mu,this%bulk,this%x_bc,this%y_bc,this%z_bc,this%dt,this%intSharp_pfloor,this%yMetric,this%dy_stretch)

        if (this%PTeqb) then
            ! subtract elastic energies to determine mixture hydrostatic energy.
            ! conductivity
            ! is assumed a function of only hydrostatic energy
            ehmix => viscwork ! use some storage space
            ehmix = this%e
            do imat = 1, this%mix%ns
                ehmix = ehmix - this%mix%material(imat)%Ys *this%mix%material(imat)%eel
            enddo
            call this%LAD%get_conductivity(this%rho,this%p,ehmix,this%T,this%sos,this%kap,this%x_bc,this%y_bc,this%z_bc,this%intSharp_tfloor)
        end if

       ! call this%get_tau( duidxj, viscwork )
        ! Now, associate the pointers to understand what's going on better
       ! tauxx => duidxj(:,:,:,tauxxidx); tauxy => duidxj(:,:,:,tauxyidx); tauxz => duidxj(:,:,:,tauxzidx);
       ! tauyy => duidxj(:,:,:,tauyyidx); tauyz => duidxj(:,:,:,tauyzidx); tauzz => duidxj(:,:,:,tauzzidx);
!print '(a,9(e21.14,1x))', 'dudx ', duidxj(89,1,1,1:9)
!print '(a,9(e21.14,1x))', 'tauxx', tauxx(89,1,1)
        ! Add the deviatoric stress to the tau for use in fluxes
        !tauxx = tauxx + this%sxx; tauxy = tauxy + this%sxy; tauxz = tauxz +this%sxz
        !                          tauyy = tauyy + this%syy; tauyz = tauyz +this%syz
        !                                                    tauzz = tauzz + this%szz
!print '(a,9(e21.14,1x))', 'tauxx', tauxx(89,1,1)

        ! store artificial stress tensor in devstress. this should not break
        ! anything since devstress will be
        ! overwritten in get_primitive and post_bc. used in update_eh -- NSG
       ! this%sxx = tauxx - this%sxx; this%sxy = tauxy - this%sxy; this%sxz = tauxz - this%sxz
       !                              this%syy = tauyy - this%syy; this%syz = tauyz - this%syz
       !                                                           this%szz = tauzz - this%szz

        ! Get heat conduction vector (q). Stored in remaining 3 components of
        ! duidxj
!        qx => duidxj(:,:,:,qxidx); qy => duidxj(:,:,:,qyidx); qz =>duidxj(:,:,:,qzidx);
!        call this%mix%get_qmix(qx, qy, qz)                     ! Get onlyspecies diffusion fluxes if PTeqb, else, everything
!        if (this%PTeqb) call this%get_q(qx, qy, qz)            ! add artificial thermal conduction fluxes



        rhs = zero

        ! regular convection and pressure terms
            ! x-derivatives, convective conservative form
           ! select case(this%mix%ns)
            !case(1)
           !     flux = this%Wcnsrv(:,:,:,mom_index)   ! mass
           !     call transpose_y_to_x(flux,xtmp1,this%decomp)
           !     call this%der%ddx(xtmp1,xtmp2,0,0)
           !     call transpose_x_to_y(xtmp2,flux,this%decomp)
           !     rhs(:,:,:,1) = rhs(:,:,:,1) - flux
           ! case default
           !     do i = 1,this%mix%ns
           !         flux = this%Wcnsrv(:,:,:,i)*this%u   ! mass
           !         call transpose_y_to_x(flux,xtmp1,this%decomp)
           !         call this%der%ddx(xtmp1,xtmp2,0,0)
           !         call transpose_x_to_y(xtmp2,flux,this%decomp)
           !         rhs(:,:,:,i) = rhs(:,:,:,i) - flux
           !     end do
           ! end select
            flux = this%Wcnsrv(:,:,:,mom_index)*this%u + this%p
            call transpose_y_to_x(flux,xtmp1,this%decomp)
            call this%der%ddx(xtmp1,xtmp2,0,0)
            call transpose_x_to_y(xtmp2,flux,this%decomp)
            rhs(:,:,:,mom_index) = rhs(:,:,:,mom_index) - flux
            flux = this%Wcnsrv(:,:,:,mom_index+1)*this%u
            call transpose_y_to_x(flux,xtmp1,this%decomp)
            call this%der%ddx(xtmp1,xtmp2,0,0)
            call transpose_x_to_y(xtmp2,flux,this%decomp)
            rhs(:,:,:,mom_index+1) = rhs(:,:,:,mom_index+1) - flux
            flux = this%Wcnsrv(:,:,:,mom_index+2)*this%u
            call transpose_y_to_x(flux,xtmp1,this%decomp)
            call this%der%ddx(xtmp1,xtmp2,0,0)
            call transpose_x_to_y(xtmp2,flux,this%decomp)
            rhs(:,:,:,mom_index+2) = rhs(:,:,:,mom_index+2) - flux
            flux = (this%Wcnsrv(:,:,:, TE_index) + this%p)*this%u
            call transpose_y_to_x(flux,xtmp1,this%decomp)
            call this%der%ddx(xtmp1,xtmp2,0,0)
            call transpose_x_to_y(xtmp2,flux,this%decomp)
            rhs(:,:,:, TE_index  ) = rhs(:,:,:, TE_index) - flux

            ! y-derivatives, convective conservative form
           ! select case(this%mix%ns)
           ! case(1)
           !     flux = this%Wcnsrv(:,:,:,mom_index+1)   ! mass
           !     call this%der%ddy(flux,tmp,0,0)
           !     rhs(:,:,:,1) = rhs(:,:,:,1) - tmp
           ! case default
           !     do i = 1,this%mix%ns
           !         flux = this%Wcnsrv(:,:,:,i)*this%v   ! mass
           !         call this%der%ddy(flux,tmp,0,0)
           !         rhs(:,:,:,i) = rhs(:,:,:,i) - tmp
           !     end do
           ! end select
            flux = this%Wcnsrv(:,:,:,mom_index)*this%v     ! x-momentum
            call this%der%ddy(flux,tmp,0,0)
            rhs(:,:,:,mom_index) = rhs(:,:,:,mom_index) - tmp
            flux = this%Wcnsrv(:,:,:,mom_index+1)*this%v + this%p  ! y-momentum
            call this%der%ddy(flux,tmp,0,0)
            rhs(:,:,:,mom_index+1) = rhs(:,:,:,mom_index+1) - tmp
            flux = this%Wcnsrv(:,:,:,mom_index+2)*this%v    ! z-momentum
            call this%der%ddy(flux,tmp,0,0)
            rhs(:,:,:,mom_index+2) = rhs(:,:,:,mom_index+2) - tmp
            flux = (this%Wcnsrv(:,:,:, TE_index) + this%p)*this%v  ! TotalEnergy
            call this%der%ddy(flux,tmp,0,0)
            rhs(:,:,:, TE_index) = rhs(:,:,:, TE_index) - tmp

            ! z-derivatives, convective conservative form
            !select case(this%mix%ns)
            !case(1)
            !    flux = this%Wcnsrv(:,:,:,mom_index+2)   ! mass
            !    call transpose_y_to_z(flux,ztmp1,this%decomp)
            !    call this%der%ddz(ztmp1,ztmp2,0,0)
            !    call transpose_z_to_y(ztmp2,flux,this%decomp)
            !    rhs(:,:,:,1) = rhs(:,:,:,1) - flux
            !case default
            !    do i = 1,this%mix%ns
            !        flux = this%Wcnsrv(:,:,:,i)*this%w   ! mass
            !        call transpose_y_to_z(flux,ztmp1,this%decomp)
            !        call this%der%ddz(ztmp1,ztmp2,0,0)
            !        call transpose_z_to_y(ztmp2,flux,this%decomp)
           !         rhs(:,:,:,i) = rhs(:,:,:,i) - flux
           !     end do
           ! end select
            flux = this%Wcnsrv(:,:,:,mom_index)*this%w     ! x-momentum
            call transpose_y_to_z(flux,ztmp1,this%decomp)
            call this%der%ddz(ztmp1,ztmp2,0,0)
            call transpose_z_to_y(ztmp2,flux,this%decomp)
            rhs(:,:,:,mom_index) = rhs(:,:,:,mom_index) - flux
            flux = this%Wcnsrv(:,:,:,mom_index+1)*this%w   ! y-momentum
            call transpose_y_to_z(flux,ztmp1,this%decomp)
            call this%der%ddz(ztmp1,ztmp2,0,0)
            call transpose_z_to_y(ztmp2,flux,this%decomp)
            rhs(:,:,:,mom_index+1) = rhs(:,:,:,mom_index+1) - flux
            flux = this%Wcnsrv(:,:,:,mom_index+2)*this%w + this%p  ! z-momentum
            call transpose_y_to_z(flux,ztmp1,this%decomp)
            call this%der%ddz(ztmp1,ztmp2,0,0)
            call transpose_z_to_y(ztmp2,flux,this%decomp)
            rhs(:,:,:,mom_index+2) = rhs(:,:,:,mom_index+2) - flux
            flux = (this%Wcnsrv(:,:,:, TE_index) + this%p)*this%w   ! Total
            call transpose_y_to_z(flux,ztmp1,this%decomp)
            call this%der%ddz(ztmp1,ztmp2,0,0)
            call transpose_z_to_y(ztmp2,flux,this%decomp)
            rhs(:,:,:, TE_index) = rhs(:,:,:, TE_index) - flux
       

        ! heat conduction (kap) terms, for energy equation
        dTdx => this%ybuf(:,:,:,1)
        dTdy => this%ybuf(:,:,:,2)
        dTdz => this%ybuf(:,:,:,3)
        d2Tdx2 => this%ybuf(:,:,:,4)
        d2Tdy2 => this%ybuf(:,:,:,5)
        d2Tdz2 => this%ybuf(:,:,:,6)
        dkapdx => NCbuff(:,:,:,1)
        dkapdy => NCbuff(:,:,:,2)
        dkapdz => NCbuff(:,:,:,3)
        call this%gradient(this%T,dTdx,dTdy,dTdz, this%x_bc, this%y_bc,this%z_bc)
        call this%secondder(this%T,d2Tdx2,d2Tdy2,d2Tdz2,[0,0],[0,0],[0,0])
        call this%gradient(this%kap,dkapdx,dkapdy,dkapdz, this%x_bc, this%y_bc,this%z_bc)
        rhs(:,:,:, TE_index) = rhs(:,:,:, TE_index) + dkapdx*dTdx+dkapdy*dTdy+dkapdz*dTdz+this%kap*(d2Tdx2+d2Tdy2+d2Tdz2)

        ! stress (mu/bulk) terms
        bambda => NCbuff(:,:,:,8)             ! used for tauxx,tauyy,tauzz
        lambda => NCbuff(:,:,:,9)
        bambda = (four/three)*this%mu + this%bulk
        lambda = this%bulk - (two/three)*this%mu
        dmudx => this%ybuf(:,:,:,1)
        dmudy => this%ybuf(:,:,:,2)
        dmudz => this%ybuf(:,:,:,3)
        dbulkdx => this%ybuf(:,:,:,4)
        dbulkdy => this%ybuf(:,:,:,5)
        dbulkdz => this%ybuf(:,:,:,6)
        call this%gradient(this%mu,dmudx,dmudy,dmudz, this%x_bc, this%y_bc,this%z_bc)
        call this%gradient(this%bulk,dbulkdx,dbulkdy,dbulkdz, this%x_bc,this%y_bc, this%z_bc)
        !tau_xx,tau_xy,tau_xz setup
        d2udx2 => NCbuff(:,:,:,1)
        d2udy2 => NCbuff(:,:,:,2)
        d2udz2 => NCbuff(:,:,:,3)
        call this%secondder(this%u,d2udx2,d2udy2,d2udz2,[0,0],[0,0],[0,0])
        !tau_xx in x-momentum and energy eq
        flux =(four/three*dmudx+dbulkdx)*dudx+(dbulkdx-two/three*dmudx)*(dvdy+dwdz)+bambda*d2udx2
        tmp = dvdy+dwdz
        call transpose_y_to_x(tmp,xtmp1,this%decomp)
        call this%der%ddx(xtmp1,xtmp2,0,0)
        call transpose_x_to_y(xtmp2,tmp,this%decomp)                   !tmp =
        flux = flux+lambda*tmp
        rhs(:,:,:, mom_index) = rhs(:,:,:, mom_index) + flux
        rhs(:,:,:, TE_index) = rhs(:,:,:, TE_index) + this%u*flux
        !tau_xy in x-momentum and energy eq
        call this%der%ddy(dvdx,tmp,0,0)        !tmp = ddy(dvdx)
        flux = dmudy*(dudy+dvdx)+this%mu*(d2udy2+tmp)
        rhs(:,:,:, mom_index) = rhs(:,:,:, mom_index) + flux
        rhs(:,:,:, TE_index) = rhs(:,:,:, TE_index) + this%u*flux
        !tau_xz in x-momentum and energy eq
        call transpose_y_to_z(dwdx,ztmp1,this%decomp)
        call this%der%ddz(ztmp1,ztmp2,0,0)
        call transpose_z_to_y(ztmp2,tmp,this%decomp)                   !tmp =
        flux = dmudz*(dudz+dwdx)+this%mu*(d2udz2+tmp)
        rhs(:,:,:, mom_index) = rhs(:,:,:, mom_index) + flux
        rhs(:,:,:, TE_index) = rhs(:,:,:, TE_index) + this%u*flux
        !tau_yy,tau_yx,tau_yz setup
        d2vdx2 => NCbuff(:,:,:,1)
        d2vdy2 => NCbuff(:,:,:,2)
        d2vdz2 => NCbuff(:,:,:,3)
        call this%secondder(this%v,d2vdx2,d2vdy2,d2vdz2,[0,0],[0,0],[0,0])
        !tau_yy in y-momentum and energy eq
        flux =(four/three*dmudy+dbulkdy)*dvdy+(dbulkdy-two/three*dmudy)*(dudx+dwdz)+bambda*d2vdy2
        call this%der%ddy(dudx+dwdz,tmp,0,0)     !tmp = ddy(dudx+dwdz)
        flux = flux+lambda*tmp
        rhs(:,:,:, mom_index+1) = rhs(:,:,:, mom_index+1) + flux
        rhs(:,:,:, TE_index) = rhs(:,:,:, TE_index) + this%v*flux
        !tau_yx in y-momentum and energy eq
        call transpose_y_to_x(dudy,xtmp1,this%decomp)
        call this%der%ddx(xtmp1,xtmp2,0,0)
        call transpose_x_to_y(xtmp2,tmp,this%decomp)                   !tmp =
        flux = dmudx*(dvdx+dudy)+this%mu*(d2vdx2+tmp)
        rhs(:,:,:, mom_index+1) = rhs(:,:,:, mom_index+1) + flux
        rhs(:,:,:, TE_index) = rhs(:,:,:, TE_index) + this%v*flux
        !tau_yz in y-momentum and energy eq
        call transpose_y_to_z(dwdy,ztmp1,this%decomp)
        call this%der%ddz(ztmp1,ztmp2,0,0)
        call transpose_z_to_y(ztmp2,tmp,this%decomp)                   !tmp =
        flux = dmudz*(dvdz+dwdy)+this%mu*(d2vdz2+tmp)
        rhs(:,:,:, mom_index+1) = rhs(:,:,:, mom_index+1) + flux
        rhs(:,:,:, TE_index) = rhs(:,:,:, TE_index) + this%v*flux
        !tau_zz,tau_zx,tau_zy setup
        d2wdx2 => NCbuff(:,:,:,1)
        d2wdy2 => NCbuff(:,:,:,2)
        d2wdz2 => NCbuff(:,:,:,3)
        call this%secondder(this%w,d2wdx2,d2wdy2,d2wdz2,[0,0],[0,0],[0,0])
        !tau_zz in z-momentum and energy eq
        flux = (four/three*dmudz+dbulkdz)*dwdz+(dbulkdz-two/three*dmudz)*(dudx+dvdy)+bambda*d2wdz2
        tmp = dudx+dvdy
        call transpose_y_to_z(tmp,ztmp1,this%decomp)
        call this%der%ddz(ztmp1,ztmp2,0,0)
        call transpose_z_to_y(ztmp2,tmp,this%decomp)                   !tmp =
        flux = flux+lambda*tmp
        rhs(:,:,:, mom_index+2) = rhs(:,:,:, mom_index+2) + flux
        rhs(:,:,:, TE_index) = rhs(:,:,:, TE_index) + this%w*flux
        !tau_zx in z-momentum and energy eq
        call transpose_y_to_x(dudz,xtmp1,this%decomp)
        call this%der%ddx(xtmp1,xtmp2,0,0)
        call transpose_x_to_y(xtmp2,tmp,this%decomp)                   !tmp =
        flux = dmudx*(dwdx+dudz)+this%mu*(d2wdx2+tmp)
        rhs(:,:,:, mom_index+2) = rhs(:,:,:, mom_index+2) + flux
        rhs(:,:,:, TE_index) = rhs(:,:,:, TE_index) + this%w*flux
        !tau_zy in z-momentum and energy eq
        call this%der%ddy(dvdz,tmp,0,0)        !tmp = ddy(dvdz)
        flux = dmudy*(dwdy+dvdz)+this%mu*(d2wdy2+tmp)
        rhs(:,:,:, mom_index+2) = rhs(:,:,:, mom_index+2) + flux
        rhs(:,:,:, TE_index) = rhs(:,:,:, TE_index) + this%w*flux
        !finish terms in energy eq
        NCbuff = duidxj
        dudx => NCbuff(:,:,:,1); dudy => NCbuff(:,:,:,2); dudz => NCbuff(:,:,:,3);
        dvdx => NCbuff(:,:,:,4); dvdy => NCbuff(:,:,:,5); dvdz => NCbuff(:,:,:,6);
        dwdx => NCbuff(:,:,:,7); dwdy => NCbuff(:,:,:,8); dwdz => NCbuff(:,:,:,9);
        call this%get_tau( duidxj, viscwork )
        tauxx => duidxj(:,:,:,tauxxidx); tauxy => duidxj(:,:,:,tauxyidx)
        tauxz => duidxj(:,:,:,tauxzidx); tauyy => duidxj(:,:,:,tauyyidx) 
        tauyz => duidxj(:,:,:,tauyzidx); tauzz => duidxj(:,:,:,tauzzidx)
        rhs(:,:,:, TE_index) = rhs(:,:,:, TE_index) + dudx*tauxx+dudy*tauxy+dudz*tauxz
        rhs(:,:,:, TE_index) = rhs(:,:,:, TE_index) + dvdx*tauxy+dvdy*tauyy+dvdz*tauyz
        rhs(:,:,:, TE_index) = rhs(:,:,:, TE_index) + dwdx*tauxz+dwdy*tauyz+dwdz*tauzz
     

         call this%mix%get_J(this%rho)

        if (this%mix%ns .GT. 3) then
            call GracefulExit("Only up to 3 species are currently supported",3214)
        end if
        ! calculating dJ/dx terms
        if (this%mix%ns .GT. 1) then
            d2Ysdx2 => NCbuff(:,:,:,1:this%mix%ns)
            d2Ysdy2 => NCbuff(:,:,:,this%mix%ns+1:2*this%mix%ns)
            d2Ysdz2 => NCbuff(:,:,:,2*this%mix%ns+1:3*this%mix%ns)
            dsumJxdx => this%ybuf(:,:,:,1); dsumJydy => this%ybuf(:,:,:,2);dsumJzdz => this%ybuf(:,:,:,3); 
            dJxdx => this%ybuf(:,:,:,4); dJydy => this%ybuf(:,:,:,5); dJzdz =>this%ybuf(:,:,:,6);
            do i = 1,this%mix%ns
                call this%secondder(this%mix%material(i)%Ys,d2Ysdx2(:,:,:,i),d2Ysdy2(:,:,:,i),d2Ysdz2(:,:,:,i),[0,0],[0,0],[0,0])
            end do
            !x-derivatives
            do i = 1,this%mix%ns
                dsumJxdx = zero
                if (this%mix%ns .GT. 2) then
                    do k = 1,this%mix%ns
                    ! this calculates the molecular diffusion flux correction
                        flux = this%Wcnsrv(:,:,:,i)*this%mix%material(k)%diff       
                        call transpose_y_to_x(flux,xtmp1,this%decomp)
                        call this%der%ddx(xtmp1,xtmp2,0,0)
                        call transpose_x_to_y(xtmp2,tmp,this%decomp)
                        dsumJxdx =dsumJxdx+tmp*dYsdx(:,:,:,k)+flux*d2Ysdx2(:,:,:,k)
                    end do
                end if
                flux = this%rho*this%mix%material(i)%diff
                call transpose_y_to_x(flux,xtmp1,this%decomp)
                call this%der%ddx(xtmp1,xtmp2,0,0)
                call transpose_x_to_y(xtmp2,tmp,this%decomp)
                dJxdx = tmp*dYsdx(:,:,:,i) + flux*d2Ysdx2(:,:,:,i) - dsumJxdx
                !species
                rhs(:,:,:,i) = rhs(:,:,:,i) + dJxdx
                !energy
                call this%mix%material(i)%get_enthalpy(tmp)    ! tmp= h_i
                rhs(:,:,:, TE_index) = rhs(:,:,:, TE_index) + tmp*dJxdx
            end do
            !y-derivatives
            do i = 1,this%mix%ns
                dsumJydy = zero
                if (this%mix%ns .GT. 2) then
                    do k = 1,this%mix%ns
                    ! this calculates the molecular diffusion flux correction
                        flux = this%Wcnsrv(:,:,:,i)*this%mix%material(k)%diff
                        call this%der%ddy(flux,tmp,0,0)
                        dsumJydy = dsumJydy+tmp*dYsdy(:,:,:,k)+flux*d2Ysdy2(:,:,:,k)
                    end do
                end if
                flux = this%rho*this%mix%material(i)%diff
                call this%der%ddy(flux,tmp,0,0)
                dJydy = tmp*dYsdy(:,:,:,i) + flux*d2Ysdy2(:,:,:,i) - dsumJydy
                !species
                rhs(:,:,:,i) = rhs(:,:,:,i) + dJydy
                !energy
                call this%mix%material(i)%get_enthalpy(tmp)    ! tmp= h_i
                rhs(:,:,:, TE_index) = rhs(:,:,:, TE_index) + tmp*dJydy
            end do
            !z-derivatives
            do i = 1,this%mix%ns
                dsumJzdz = zero
                if (this%mix%ns .GT. 2) then
                    do k = 1,this%mix%ns
                    ! this calculates the molecular diffusion flux correction
                        flux = this%Wcnsrv(:,:,:,i)*this%mix%material(k)%diff
                        call transpose_y_to_z(flux,ztmp1,this%decomp)
                        call this%der%ddz(ztmp1,ztmp2,0,0)
                        call transpose_z_to_y(ztmp2,tmp,this%decomp)
                        dsumJzdz = dsumJzdz+tmp*dYsdz(:,:,:,k)+flux*d2Ysdz2(:,:,:,k)
                    end do
                end if
                flux = this%rho*this%mix%material(i)%diff
                call transpose_y_to_z(flux,ztmp1,this%decomp)
                call this%der%ddz(ztmp1,ztmp2,0,0)
                call transpose_z_to_y(ztmp2,tmp,this%decomp)
                dJzdz = tmp*dYsdz(:,:,:,i) + flux*d2Ysdz2(:,:,:,i) - dsumJzdz
                !species
                rhs(:,:,:,i) = rhs(:,:,:,i) + dJzdz
                !energy
                call this%mix%material(i)%get_enthalpy(tmp)    ! tmp= h_i
                rhs(:,:,:, TE_index) = rhs(:,:,:, TE_index) + tmp*dJzdz
            end do

            !final J terms in energy equation:
          ! call this%get_J(gradYs)
          ! Jx => gradYs(:,:,:,1:this%mix%ns)
          ! Jy => gradYs(:,:,:,this%mix%ns+1:2*this%mix%ns)
          ! Jz => gradYs(:,:,:,2*this%mix%ns+1:3*this%mix%ns)
            do i = 1,this%mix%ns
                call this%mix%material(i)%get_enthalpy(flux)
                call transpose_y_to_x(flux,xtmp1,this%decomp)
                call this%der%ddx(xtmp1,xtmp2,this%x_bc(1),this%x_bc(2))
                call transpose_x_to_y(xtmp2,tmp,this%decomp)
                rhs(:,:,:, TE_index) = rhs(:,:,:, TE_index) +tmp*this%mix%material(i)%Ji(:,:,:,1)
                call this%der%ddy(flux,tmp,this%y_bc(1),this%y_bc(2))
                rhs(:,:,:, TE_index) = rhs(:,:,:, TE_index) + tmp*this%mix%material(i)%Ji(:,:,:,2)
                call transpose_y_to_z(flux,ztmp1,this%decomp)
                call this%der%ddz(ztmp1,ztmp2,this%z_bc(1),this%z_bc(2))
                call transpose_z_to_y(ztmp2,tmp,this%decomp)
                rhs(:,:,:, TE_index) = rhs(:,:,:, TE_index) + tmp*this%mix%material(i)%Ji(:,:,:,3)
            end do
        end if 


       if(this%intSharp.AND.this%intSharp_cpl) then
           !calculate kinetic energy for intSharp terms
           ke = half*( this%u**2 + this%v**2 + this%w**2 ) !is this accesiblewithout recreating?

           if(this%intSharp_spf) then
              ! if(useNewSPF) then !this is unstable
              !    !new -- for useNewSPF = .TRUE. in solidmix
                 rhs(:,:,:,mom_index  ) = rhs(:,:,:,mom_index  ) + this%mix%intSharp_f(:,:,:,1)
                 rhs(:,:,:,mom_index+1) = rhs(:,:,:,mom_index+1) + this%mix%intSharp_f(:,:,:,2)
                 rhs(:,:,:,mom_index+2) = rhs(:,:,:,mom_index+2) + this%mix%intSharp_f(:,:,:,3)
                 rhs(:,:,:,TE_index   ) = rhs(:,:,:,TE_index   ) + this%mix%intSharp_h(:,:,:,1)
              ! else
              !    !original
              !    rhs(:,:,:,mom_index  ) = rhs(:,:,:,mom_index  ) +
              !    this%mix%intSharp_f(:,:,:,1)*this%u
              !    rhs(:,:,:,mom_index+1) = rhs(:,:,:,mom_index+1) +
              !    this%mix%intSharp_f(:,:,:,1)*this%v
              !    rhs(:,:,:,mom_index+2) = rhs(:,:,:,mom_index+2) +
              !    this%mix%intSharp_f(:,:,:,1)*this%w
              !    rhs(:,:,:,TE_index   ) = rhs(:,:,:,TE_index   ) +
              !    this%mix%intSharp_f(:,:,:,1)*ke +
              !    this%mix%intSharp_h(:,:,:,1)
              ! endif

              !high order VF bounds diffusion terms
              call divergence(this%decomp,this%der,this%mix%intSharp_fDiff(:,:,:,1)*this%u,this%mix%intSharp_fDiff(:,:,:,2)*this%u,this%mix%intSharp_fDiff(:,:,:,3)*this%u,tmp,this%x_bc,-this%y_bc,-this%z_bc)
              rhs(:,:,:,mom_index  ) = rhs(:,:,:,mom_index  ) + tmp

              call divergence(this%decomp,this%der,this%mix%intSharp_fDiff(:,:,:,1)*this%v,this%mix%intSharp_fDiff(:,:,:,2)*this%v,this%mix%intSharp_fDiff(:,:,:,3)*this%v,tmp,-this%x_bc,this%y_bc,-this%z_bc)
              rhs(:,:,:,mom_index+1) = rhs(:,:,:,mom_index+1) + tmp

              call divergence(this%decomp,this%der,this%mix%intSharp_fDiff(:,:,:,1)*this%w,this%mix%intSharp_fDiff(:,:,:,2)*this%w,this%mix%intSharp_fDiff(:,:,:,3)*this%w,tmp,-this%x_bc,-this%y_bc,this%z_bc)
              rhs(:,:,:,mom_index+2) = rhs(:,:,:,mom_index+2) + tmp

              call divergence(this%decomp,this%der,this%mix%intSharp_fDiff(:,:,:,1)*ke + this%mix%intSharp_hDiff(:,:,:,1),this%mix%intSharp_fDiff(:,:,:,2)*ke + this%mix%intSharp_hDiff(:,:,:,2),this%mix%intSharp_fDiff(:,:,:,3)*ke + this%mix%intSharp_hDiff(:,:,:,3),tmp,-this%x_bc,-this%y_bc,-this%z_bc)
              rhs(:,:,:,TE_index   ) = rhs(:,:,:,TE_index   ) + tmp

           else

              !low order terms
              call divergence(this%decomp,this%derD02,this%mix%intSharp_f(:,:,:,1)*this%u,this%mix%intSharp_f(:,:,:,2)*this%u,this%mix%intSharp_f(:,:,:,3)*this%u,tmp,this%x_bc,-this%y_bc,-this%z_bc)
              rhs(:,:,:,mom_index  ) = rhs(:,:,:,mom_index  ) + tmp

              call divergence(this%decomp,this%derD02,this%mix%intSharp_f(:,:,:,1)*this%v,this%mix%intSharp_f(:,:,:,2)*this%v,this%mix%intSharp_f(:,:,:,3)*this%v,tmp,-this%x_bc,this%y_bc,-this%z_bc)
              rhs(:,:,:,mom_index+1) = rhs(:,:,:,mom_index+1) + tmp

              call divergence(this%decomp,this%derD02,this%mix%intSharp_f(:,:,:,1)*this%w,this%mix%intSharp_f(:,:,:,2)*this%w,this%mix%intSharp_f(:,:,:,3)*this%w,tmp,-this%x_bc,-this%y_bc,this%z_bc)
              rhs(:,:,:,mom_index+2) = rhs(:,:,:,mom_index+2) + tmp

              call divergence(this%decomp,this%derD02,this%mix%intSharp_f(:,:,:,1)*ke + this%mix%intSharp_h(:,:,:,1),this%mix%intSharp_f(:,:,:,2)*ke + this%mix%intSharp_h(:,:,:,2),this%mix%intSharp_f(:,:,:,3)*ke + this%mix%intSharp_h(:,:,:,3),tmp,-this%x_bc,-this%y_bc,-this%z_bc)
              rhs(:,:,:,TE_index   ) = rhs(:,:,:,TE_index   ) + tmp


              !high order terms
              call divergence(this%decomp,this%der,this%mix%intSharp_fDiff(:,:,:,1)*this%u,this%mix%intSharp_fDiff(:,:,:,2)*this%u,this%mix%intSharp_fDiff(:,:,:,3)*this%u,tmp,this%x_bc,-this%y_bc,-this%z_bc)
              rhs(:,:,:,mom_index  ) = rhs(:,:,:,mom_index  ) + tmp

              call divergence(this%decomp,this%der,this%mix%intSharp_fDiff(:,:,:,1)*this%v,this%mix%intSharp_fDiff(:,:,:,2)*this%v,this%mix%intSharp_fDiff(:,:,:,3)*this%v,tmp,-this%x_bc,this%y_bc,-this%z_bc)
              rhs(:,:,:,mom_index+1) = rhs(:,:,:,mom_index+1) + tmp

              call divergence(this%decomp,this%der,this%mix%intSharp_fDiff(:,:,:,1)*this%w,this%mix%intSharp_fDiff(:,:,:,2)*this%w,this%mix%intSharp_fDiff(:,:,:,3)*this%w,tmp,-this%x_bc,-this%y_bc,this%z_bc)
              rhs(:,:,:,mom_index+2) = rhs(:,:,:,mom_index+2) + tmp

              call divergence(this%decomp,this%der,this%mix%intSharp_fDiff(:,:,:,1)*ke + this%mix%intSharp_hDiff(:,:,:,1),this%mix%intSharp_fDiff(:,:,:,2)*ke + this%mix%intSharp_hDiff(:,:,:,2),this%mix%intSharp_fDiff(:,:,:,3)*ke + this%mix%intSharp_hDiff(:,:,:,3),tmp,-this%x_bc,-this%y_bc,-this%z_bc)
              rhs(:,:,:,TE_index   ) = rhs(:,:,:,TE_index   ) + tmp

              !FV sharpening
              rhs(:,:,:,mom_index  ) = rhs(:,:,:,mom_index  ) + this%mix%intSharp_fFV(:,:,:,1)
              rhs(:,:,:,mom_index+1) = rhs(:,:,:,mom_index+1) + this%mix%intSharp_fFV(:,:,:,2)
              rhs(:,:,:,mom_index+2) = rhs(:,:,:,mom_index+2) + this%mix%intSharp_fFV(:,:,:,3)
              rhs(:,:,:,TE_index   ) = rhs(:,:,:,TE_index   ) + this%mix%intSharp_hFV

           endif
    endif
        


        !!Surface Tension
        if (this%use_surfaceTension) then
            rhs(:,:,:,mom_index  ) = rhs(:,:,:,mom_index  ) + this%mix%surfaceTension_f(:,:,:,1)
            rhs(:,:,:,mom_index+1) = rhs(:,:,:,mom_index+1) + this%mix%surfaceTension_f(:,:,:,2)
            rhs(:,:,:,mom_index+2) = rhs(:,:,:,mom_index+2) + this%mix%surfaceTension_f(:,:,:,3)
            rhs(:,:,:,TE_index   ) = rhs(:,:,:,TE_index   ) + this%mix%surfaceTension_e
        endif

        ! Call problem source hook
        
       ! Call problem source hook
       call hook_mixture_source(this%decomp, this%mesh, this%fields, this%mix, this%tsim, rhs)




        call this%get_tau( duidxj, viscwork )
        ! Now, associate the pointers to understand what's going on better
        tauxx => duidxj(:,:,:,tauxxidx); tauxy => duidxj(:,:,:,tauxyidx); tauxz => duidxj(:,:,:,tauxzidx);
        tauyy => duidxj(:,:,:,tauyyidx); tauyz => duidxj(:,:,:,tauyzidx); tauzz => duidxj(:,:,:,tauzzidx);
!print '(a,9(e21.14,1x))', 'dudx ', duidxj(89,1,1,1:9)
!print '(a,9(e21.14,1x))', 'tauxx', tauxx(89,1,1)
        ! Add the deviatoric stress to the tau for use in fluxes
        tauxx = tauxx + this%sxx; tauxy = tauxy + this%sxy; tauxz = tauxz + this%sxz
                                  tauyy = tauyy + this%syy; tauyz = tauyz + this%syz
                                                            tauzz = tauzz + this%szz
!print '(a,9(e21.14,1x))', 'tauxx', tauxx(89,1,1)

        ! store artificial stress tensor in devstress. this should not break
        ! anything since devstress will be
        ! overwritten in get_primitive and post_bc. used in update_eh -- NSG
        this%sxx = tauxx - this%sxx; this%sxy = tauxy - this%sxy; this%sxz =tauxz - this%sxz
                                     this%syy = tauyy - this%syy; this%syz =tauyz - this%syz
                                                                  this%szz =tauzz - this%szz

    end subroutine

    subroutine getRHS_xStagg( this,  rhs, tauxx,tauxy,tauxz, qx)
        use operators, only: gradFV_x, interpolateFV_x, gradFV_N2Fx
        class(sgrid), target, intent(inout) :: this
        real(rkind), dimension(this%nxp, this%nyp, this%nzp, ncnsrv),intent(inout) :: rhs
        real(rkind), dimension(this%nxp, this%nyp, this%nzp), intent(in) :: tauxx,tauxy,tauxz
        real(rkind), dimension(this%nxp, this%nyp, this%nzp), intent(in) :: qx
        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: buff, flux,TE, u_int,v_int, w_int, p_int, tauxx_int, tauxy_int, tauxz_int, qx_int, e_int, rho_int, rhodiff_int
        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: rhou_int,rhov_int, rhow_int, rhoe_int,spe_int, rhoYs_int, den, num, gradRYs, tauRho_mid, rhom, rhom_int
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,2) :: VF_int 
        real(rkind), dimension(:,:,:), pointer :: xtmp1,xtmp2
        integer :: i

        call interpolateFV_x(this%decomp,this%interpMid,this%u,u_int,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        call interpolateFV_x(this%decomp,this%interpMid,this%v,v_int,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        call interpolateFV_x(this%decomp,this%interpMid,this%w,w_int,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        call interpolateFV_x(this%decomp,this%interpMid,this%p,p_int,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        call interpolateFV_x(this%decomp,this%interpMid,this%rho,rho_int,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        call interpolateFV_x(this%decomp,this%interpMid,this%e,e_int,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        call interpolateFV_x(this%decomp,this%interpMid,qx,qx_int,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)

        call interpolateFV_x(this%decomp,this%interpMid,this%rho*this%e,rhoe_int,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)

        tauRho_mid = 0.0
        do i = 1,2

          call this%mix%material(i)%getSpeciesDensity(this%rho,rhom)
          call interpolateFV_x(this%decomp,this%interpMid,rhom,rhom_int,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
          call interpolateFV_x(this%decomp,this%interpMid,this%mix%material(i)%rhodiff,rhodiff_int,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
          call gradFV_N2Fx(this%decomp,this%derStagg,this%mix%material(i)%Ys*this%rho,gradRYs,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)

          tauRho_mid = tauRho_mid + tauxx/rho_int*rhodiff_int*gradRYs


        enddo

        call gradFV_x(this%decomp,this%derStagg,tauRho_mid,this%tauRho,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        !rho_int = 0.0
        !num = 0.0; den = 0.0;
!        call interpolateFV_x(this%decomp,this%interpMid,this%mix%material(1)%VF,VF_int,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
!        VF_int(:,:,:,2) = 1 - VF_int(:,:,:,1)

        !do i = 1,2

        !  call interpolateFV_x(this%decomp,this%interpMid,this%rho*this%mix%material(i)%Ys,rhoYs_int,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        !  call interpolateFV_x(this%decomp,this%interpMid,this%mix%material(i)%VF,VF_int(:,:,:,i),this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        !  rho_int = rho_int + rhoYs_int
        !  num = num + VF_int(:,:,:,i)*this%mix%material(i)%hydro%gam*this%mix%material(i)%hydro%Pinf*this%mix%material(i)%hydro%onebygam_m1
        !  den = den + VF_int(:,:,:,i)*this%mix%material(i)%hydro%onebygam_m1
        !enddo

        !e_int = (1/rho_int)*(p_int*den + num)

        flux = 0.0
        buff = rho_int*u_int*u_int + p_int - tauxx !x-momentum
        !buff = rhou_int*u_int + p_int - tauxx_int
!print *, 'flux 1', flux(89,1,1), this%u(89,1,1), this%p(89,1,1), tauxx(89,1,1)
        !if(this%use_CnsrvSurfaceTension) then

        !    buff = buff - this%mix%surfaceTension_fxx

        !endif
        call gradFV_x(this%decomp,this%derStagg,buff,flux,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        rhs(:,:,:,mom_index  ) = rhs(:,:,:,mom_index  ) - flux
        this%xflux_x = flux

        flux =0.0
        buff = rho_int*u_int*v_int  - tauxy !y-momentum
        !buff = rhov_int*u_int - tauxy_int
        !if(this%use_CnsrvSurfaceTension) then

        !    buff = buff - this%mix%surfaceTension_fxy_x

        !endif
        call gradFV_x(this%decomp,this%derStagg,buff,flux,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        rhs(:,:,:,mom_index+1) = rhs(:,:,:,mom_index+1) - flux
        this%xflux_y = flux

        buff = rho_int*u_int*w_int - tauxz !z-momentum
        !if(this%use_CnsrvSurfaceTension) then

        !    buff = buff - this%mix%surfaceTension_fxz_x

        !endif
        !buff = rhow_int*u_int - tauxz_int
        call gradFV_x(this%decomp,this%derStagg,buff,flux,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        rhs(:,:,:,mom_index+2) = rhs(:,:,:,mom_index+2) - flux
        this%xflux_z = flux
        flux = 0.0

        TE = rhoe_int +rho_int*half*(u_int*u_int + v_int*v_int + w_int*w_int) 
        !TE = rhoe_int + half*rho_int*(u_int*u_int + v_int*v_int + w_int*w_int)
        buff = ( TE + p_int - tauxx )*u_int - v_int*tauxy - w_int*tauxz + qx !+ tauRho_mid
        !buff = ( TE + p_int )*u_int ! qx_int - v_int*tauxy_int - w_int*tauxz_int

        !if(this%use_CnsrvSurfaceTension) then

        !    buff = buff + this%mix%surfaceTension_pe_x*u_int -this%mix%surfaceTension_fxx*u_int - this%mix%surfaceTension_fxy_x*v_int - this%mix%surfaceTension_fxz_x*w_int

        !endif       
        call gradFV_x(this%decomp,this%derStagg,buff,flux,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        rhs(:,:,:, TE_index  ) = rhs(:,:,:, TE_index  ) - flux
        this%xflux_e = flux

        this%rho_int = rho_int
        this%u_int = u_int
        this%v_int = v_int
        this%w_int = w_int
        this%e_int     = e_int
        this%p_int     = p_int

    end subroutine

    subroutine getRHS_yStagg( this,  rhs, tauxy,tauyy,tauyz, qy)
        use operators, only: gradFV_y, interpolateFV_y
        class(sgrid), target, intent(inout) :: this
        real(rkind), dimension(this%nxp, this%nyp, this%nzp,ncnsrv),intent(inout) :: rhs
        real(rkind), dimension(this%nxp, this%nyp, this%nzp), intent(in) :: tauxy,tauyy,tauyz
        real(rkind), dimension(this%nxp, this%nyp, this%nzp), intent(in) :: qy
        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: buff, flux,TE,u_int,v_int, w_int, p_int, tauxy_int, tauyy_int, tauyz_int, qy_int, e_int,rho_int
        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: rhou_int,rhov_int,rhow_int, rhoe_int, spe_int, rhoYs_int, VF_int, den, num
        real(rkind), dimension(:,:,:), pointer :: xtmp1,xtmp2
        integer :: i

        call interpolateFV_y(this%decomp,this%interpMid,this%u,u_int,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        call interpolateFV_y(this%decomp,this%interpMid,this%v,v_int,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        call interpolateFV_y(this%decomp,this%interpMid,this%w,w_int,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        call interpolateFV_y(this%decomp,this%interpMid,this%p,p_int,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        call interpolateFV_y(this%decomp,this%interpMid,this%rho,rho_int,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        call interpolateFV_y(this%decomp,this%interpMid,this%e,e_int,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        call interpolateFV_y(this%decomp,this%interpMid,qy,qy_int,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        call interpolateFV_y(this%decomp,this%interpMid,this%rho*this%e,rhoe_int,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)

        !rho_int = 0.0
        !num = 0.0; den = 0.0;
        !do i = 1,2

        !  call interpolateFV_y(this%decomp,this%interpMid,this%rho*this%mix%material(i)%Ys,rhoYs_int,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        !  call interpolateFV_y(this%decomp,this%interpMid,this%mix%material(i)%VF,VF_int,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        !  rho_int = rho_int + rhoYs_int
        !  num = num + VF_int*this%mix%material(i)%hydro%gam*this%mix%material(i)%hydro%Pinf*this%mix%material(i)%hydro%onebygam_m1
        !  den = den + VF_int*this%mix%material(i)%hydro%onebygam_m1
        !enddo
  
        !e_int = (1/rho_int)*(p_int*den + num)

        flux = 0.0
        buff = rho_int*v_int*u_int   - tauxy !x-momentum
        !buff = rhou_int*v_int - tauxy_int
!print *, 'flux 1', flux(89,1,1), this%u(89,1,1), this%p(89,1,1), tauxx(89,1,1)
        !if(this%use_CnsrvSurfaceTension) then

        !    buff = buff - this%mix%surfaceTension_fxy_y

        !endif
        call gradFV_y(this%decomp,this%derStagg,buff,flux,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        rhs(:,:,:,mom_index  ) = rhs(:,:,:,mom_index  ) - flux
        this%yflux_x = flux

        flux = 0.0
        buff = rho_int*v_int*v_int + p_int   - tauyy !y-momentum
        !buff = rhov_int*v_int + p_int - tauyy_int
        !if(this%use_CnsrvSurfaceTension) then

        !    buff = buff - this%mix%surfaceTension_fyy

        !endif
        call gradFV_y(this%decomp,this%derStagg,buff,flux,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        rhs(:,:,:,mom_index+1) = rhs(:,:,:,mom_index+1)  - flux
        this%yflux_y = flux

        flux = 0.0
        buff = rho_int*v_int*w_int   - tauyz !z-momentum
        !buff = rhow_int*w_int -tauyz_int
        !        if(this%use_CnsrvSurfaceTension) then

        !    buff = buff - this%mix%surfaceTension_fxy_y

        !endif
        call gradFV_y(this%decomp,this%derStagg,buff,flux,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        rhs(:,:,:,mom_index+2) = rhs(:,:,:,mom_index+2) - flux
        this%yflux_z = flux

        flux = 0.0
        TE = rhoe_int + rho_int*half*(u_int*u_int + v_int*v_int + w_int*w_int) 
        !TE = rhoe_int + half*(rhou_int*u_int + v_int*rhov_int + w_int*rhow_int)
        buff = ( TE + p_int - tauyy )*v_int  - u_int*tauxy -w_int*tauyz + qy
        !buff = ( TE + p_int  )*v_int 

        ! if(this%use_CnsrvSurfaceTension) then

        !    buff = buff + this%mix%surfaceTension_pe_y*v_int - this%mix%surfaceTension_fxy_y*u_int - this%mix%surfaceTension_fyy*v_int - this%mix%surfaceTension_fyz_y*w_int

        !endif

        call gradFV_y(this%decomp,this%derStagg,buff,flux,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        rhs(:,:,:, TE_index  ) = rhs(:,:,:, TE_index  ) - flux
        this%yflux_e = flux

        !this%rho_int = rho_int
        !this%u_int = u_int
        !this%v_int = v_int
        !this%w_int = w_int
        !this%e_int     = e_int
        !this%p_int     = p_int

    end subroutine

   
    subroutine getRHS_zStagg( this,  rhs, tauxz,tauyz,tauzz, qz)
        use operators, only: gradFV_z, interpolateFV_z
        class(sgrid), target, intent(inout) :: this
        real(rkind), dimension(this%nxp, this%nyp,this%nzp,ncnsrv),intent(inout) :: rhs
        real(rkind), dimension(this%nxp, this%nyp, this%nzp), intent(in) ::tauxz,tauyz,tauzz
        real(rkind), dimension(this%nxp, this%nyp, this%nzp), intent(in) :: qz
        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: buff,flux,TE,u_int,v_int, w_int, p_int, tauxz_int, tauzz_int, tauyz_int, qz_int,e_int,rho_int
        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: rhou_int,rhov_int, rhow_int, rhoe_int, spe_int, rhoYs_int, VF_int, den, num
        real(rkind), dimension(:,:,:), pointer :: xtmp1,xtmp2
        integer :: i

        call interpolateFV_z(this%decomp,this%interpMid,this%u,u_int,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        call interpolateFV_z(this%decomp,this%interpMid,this%v,v_int,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        call interpolateFV_z(this%decomp,this%interpMid,this%w,w_int,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        call interpolateFV_z(this%decomp,this%interpMid,this%p,p_int,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        call interpolateFV_z(this%decomp,this%interpMid,this%rho,rho_int,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        call interpolateFV_z(this%decomp,this%interpMid,this%e,e_int,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        call interpolateFV_z(this%decomp,this%interpMid,this%mix%surfaceTension_pe,spe_int,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        call interpolateFV_z(this%decomp,this%interpMid,qz,qz_int,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        
        !rho_int = 0.0
        !num = 0.0; den = 0.0;
        !do i = 1,2

        !  call interpolateFV_z(this%decomp,this%interpMid,this%rho*this%mix%material(i)%Ys,rhoYs_int,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        !  call interpolateFV_z(this%decomp,this%interpMid,this%mix%material(i)%VF,VF_int,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        !  rho_int = rho_int + rhoYs_int
        !  num = num + VF_int*this%mix%material(i)%hydro%gam*this%mix%material(i)%hydro%Pinf*this%mix%material(i)%hydro%onebygam_m1
        !  den = den + VF_int*this%mix%material(i)%hydro%onebygam_m1
        !enddo

        !e_int = (1/rho_int)*(p_int*den + num)

        flux  = 0
        buff  = 0
        buff = rho_int*w_int*u_int  - tauxz !x-momentum
        !buff = rhou_int*w_int -tauxz_int
        if(this%use_CnsrvSurfaceTension) then

            buff = buff - this%mix%surfaceTension_fxz_z

        endif
!print *, 'flux 1', flux(89,1,1), this%u(89,1,1), this%p(89,1,1), tauxx(89,1,1)
        call gradFV_z(this%decomp,this%derStagg,buff,flux,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        rhs(:,:,:,mom_index  ) = rhs(:,:,:,mom_index  ) - flux
        this%zflux_x = flux 

        flux = 0.0
        buff = rho_int*w_int*v_int   - tauyz !y-momentum
        !buff = rhov_int*w_int - tauyz_int
        if(this%use_CnsrvSurfaceTension) then

            buff = buff - this%mix%surfaceTension_fyz_z

        endif
        call gradFV_z(this%decomp,this%derStagg,buff,flux,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        rhs(:,:,:,mom_index+1) = rhs(:,:,:,mom_index+1) - flux
        this%zflux_y = flux

        flux = 0.0
        buff = rho_int*w_int*w_int + p_int   - tauzz !z-momentum
        !buff = rhow_int*w_int + p_int -tauzz_int
        if(this%use_CnsrvSurfaceTension) then

            buff = buff - this%mix%surfaceTension_fzz

        endif

        call gradFV_z(this%decomp,this%derStagg,buff,flux,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        rhs(:,:,:,mom_index+2) = rhs(:,:,:,mom_index+2) - flux
        this%zflux_z = flux
        flux = 0.0
        !TE = rhoe_int + half*(rhou_int*u_int + rhov_int*v_int + rhow_int*w_int)
        TE = rho_int*(e_int + half*(u_int*u_int + v_int*v_int + w_int*w_int) )

        buff = ( TE + p_int - tauzz )*w_int - u_int*tauxz-v_int*tauyz + qz
        !buff = ( TE + p_int )*w_int 
        if(this%use_CnsrvSurfaceTension) then

            buff = buff + this%mix%surfaceTension_pe_z*w_int - this%mix%surfaceTension_fxz_z*u_int - this%mix%surfaceTension_fyz_z*v_int - this%mix%surfaceTension_fzz*w_int

        endif
        call gradFV_z(this%decomp,this%derStagg,buff,flux,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        rhs(:,:,:, TE_index  ) = rhs(:,:,:, TE_index  ) - flux
        this%zflux_e = flux


    end subroutine
    
    subroutine getRHS_x( this,  rhs, tauxx,tauxy,tauxz, qx)
        use operators, only: divergence,divergenceFV,interpolateFV,interpolateFV_x,interpolateFV_y,interpolateFV_z,gradFV_y, gradFV_z, gradFV_x
        class(sgrid), target, intent(inout) :: this
        real(rkind), dimension(this%nxp, this%nyp, this%nzp, ncnsrv), intent(inout) :: rhs
        real(rkind), dimension(this%nxp, this%nyp, this%nzp), intent(in) :: tauxx,tauxy,tauxz
        real(rkind), dimension(this%nxp, this%nyp, this%nzp), intent(in) :: qx
        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: flux, ddx_p,p_int, u_int, ddx_up, ke
        real(rkind), dimension(:,:,:), pointer :: xtmp1,xtmp2
        integer :: i

        xtmp1 => this%xbuf(:,:,:,1); xtmp2 => this%xbuf(:,:,:,2)

        !call interpolateFV_x(this%decomp,this%interpMid,this%p,p_int,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        !call gradFV_x(this%decomp,this%derStagg,p_int,ddx_p,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        !flux = this%Wcnsrv(:,:,:,mom_index  )*this%u - tauxx
      
        flux = this%Wcnsrv(:,:,:,mom_index  )*this%u + this%p - tauxx  ! x-momentum
!print *, 'flux 1', flux(89,1,1), this%u(89,1,1), this%p(89,1,1), tauxx(89,1,1)
        if(this%use_CnsrvSurfaceTension) then

            flux = flux - this%mix%surfaceTension_fxx

        endif
        call transpose_y_to_x(flux,xtmp1,this%decomp)
        call this%der%ddx(xtmp1,xtmp2, this%x_bc(1), this%x_bc(2)) ! Symmetric for x-momentum
!do i = 1, size(flux,1)
!  write(*,'(4(e21.14,1x))') xtmp1(i,1,1), xtmp2(i,1,1)
!enddo
        call transpose_x_to_y(xtmp2,flux,this%decomp)
        this%xflux_x = flux

        rhs(:,:,:,mom_index  ) = rhs(:,:,:,mom_index  ) - flux 
        flux = this%Wcnsrv(:,:,:,mom_index  )*this%v          - tauxy   ! y-momentum
        if(this%use_CnsrvSurfaceTension) then

            flux = flux - this%mix%surfaceTension_fxy

        endif

        call transpose_y_to_x(flux,xtmp1,this%decomp)
        call this%der%ddx(xtmp1,xtmp2,-this%x_bc(1),-this%x_bc(2)) ! Anti-symmetric for all but x-momentum
        call transpose_x_to_y(xtmp2,flux,this%decomp)
        rhs(:,:,:,mom_index+1) = rhs(:,:,:,mom_index+1) - flux
        this%xflux_y = flux 

        flux = this%Wcnsrv(:,:,:,mom_index  )*this%w          - tauxz  ! z-momentum
        if(this%use_CnsrvSurfaceTension) then

            flux = flux - this%mix%surfaceTension_fxz

        endif

        call transpose_y_to_x(flux,xtmp1,this%decomp)
        call this%der%ddx(xtmp1,xtmp2,-this%x_bc(1),-this%x_bc(2)) ! Anti-symmetric for all but x-momentum
        call transpose_x_to_y(xtmp2,flux,this%decomp)
        rhs(:,:,:,mom_index+2) = rhs(:,:,:,mom_index+2) - flux
        this%xflux_z = flux

        !call interpolateFV_x(this%decomp,this%interpMid,this%u,u_int,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        !call gradFV_x(this%decomp,this%derStagg,p_int*u_int,ddx_up,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        !flux = (this%Wcnsrv(:,:,:, TE_index  ) - tauxx)*this%u -this%v*tauxy - this%w*tauxz + qx
        ke = 0.5*(this%u*this%u + this%v*this%v + this%w*this%w)
        flux = (this%Wcnsrv(:,:,:, TE_index  ) + this%p - tauxx)*this%u - this%v*tauxy - this%w*tauxz + qx ! Total Energy
        if(this%use_CnsrvSurfaceTension) then

            flux = flux - this%mix%surfaceTension_fxx*this%u -this%mix%surfaceTension_fxy*this%v - this%mix%surfaceTension_fxz*this%w

        endif
        call transpose_y_to_x(flux,xtmp1,this%decomp)
        call this%der%ddx(xtmp1,xtmp2,-this%x_bc(1),-this%x_bc(2)) ! Anti-symmetric for all but x-momentum
        call transpose_x_to_y(xtmp2,flux,this%decomp)
        rhs(:,:,:, TE_index  ) = rhs(:,:,:, TE_index  ) - flux 
        this%xflux_e = flux

    end subroutine

    subroutine getRHS_y( this,  rhs,&
                        tauxy,tauyy,tauyz,&
                            qy )
        use operators, only: divergence,divergenceFV,interpolateFV,interpolateFV_x,interpolateFV_y,interpolateFV_z,gradFV_y, gradFV_z, gradFV_x
        class(sgrid), target, intent(inout) :: this
        real(rkind), dimension(this%nxp, this%nyp, this%nzp, ncnsrv), intent(inout) :: rhs
        real(rkind), dimension(this%nxp, this%nyp, this%nzp), intent(in) :: tauxy,tauyy,tauyz
        real(rkind), dimension(this%nxp, this%nyp, this%nzp), intent(in) :: qy

        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: flux, p_int,ddy_p, v_int, ddy_vp
        real(rkind), dimension(:,:,:), pointer :: ytmp1

        ytmp1 => this%ybuf(:,:,:,6)

        flux = this%Wcnsrv(:,:,:,mom_index+1)*this%u          - tauxy ! x-momentum
        if(this%use_CnsrvSurfaceTension) then

            flux = flux - this%mix%surfaceTension_fxy

        endif

        call this%der%ddy(flux,ytmp1,-this%y_bc(1),-this%y_bc(2)) ! Anti-symmetric for all but y-momentum
        rhs(:,:,:,mom_index  ) = rhs(:,:,:,mom_index  ) - ytmp1
        this%yflux_x = flux

       ! call interpolateFV_y(this%decomp,this%interpMid,this%p,p_int,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
       ! call gradFV_y(this%decomp,this%derStagg,p_int,ddy_p,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        flux = this%Wcnsrv(:,:,:,mom_index+1)*this%v + this%p - tauyy ! y-momentum
        !flux = this%Wcnsrv(:,:,:,mom_index+1)*this%v  - tauyy
        if(this%use_CnsrvSurfaceTension) then

            flux = flux - this%mix%surfaceTension_fyy

        endif
        call this%der%ddy(flux,ytmp1, this%y_bc(1), this%y_bc(2)) ! Symmetric for y-momentum
        rhs(:,:,:,mom_index+1) = rhs(:,:,:,mom_index+1) - ytmp1 
        this%yflux_y = flux

        flux = this%Wcnsrv(:,:,:,mom_index+1)*this%w          - tauyz !z-momentum
        if(this%use_CnsrvSurfaceTension) then

            flux = flux - this%mix%surfaceTension_fyz

        endif
        call this%der%ddy(flux,ytmp1,-this%y_bc(1),-this%y_bc(2)) ! Anti-symmetric for all but y-momentum
        rhs(:,:,:,mom_index+2) = rhs(:,:,:,mom_index+2) - ytmp1 
        this%yflux_z = flux

        !call interpolateFV_y(this%decomp,this%interpMid,this%v,v_int,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        !call gradFV_y(this%decomp,this%derStagg,p_int*v_int,ddy_vp,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        !flux = (this%Wcnsrv(:,:,:, TE_index  ) - tauyy)*this%v - this%u*tauxy -this%w*tauyz + qy 
        flux = (this%Wcnsrv(:,:,:, TE_index  ) + this%p - tauyy)*this%v - this%u*tauxy - this%w*tauyz + qy ! Total Energy
        if(this%use_CnsrvSurfaceTension) then

            flux = flux - this%mix%surfaceTension_fxy*this%u - this%mix%surfaceTension_fyy*this%v - this%mix%surfaceTension_fyz*this%w  

        endif
        
        call this%der%ddy(flux,ytmp1,-this%y_bc(1),-this%y_bc(2)) ! Anti-symmetric for all but y-momentum
        rhs(:,:,:, TE_index  ) = rhs(:,:,:, TE_index  ) - ytmp1 
        this%yflux_e = flux

    end subroutine

    subroutine getRHS_z(       this,  rhs,&
                        tauxz,tauyz,tauzz,&
                            qz )
        use operators, only: divergence,divergenceFV,interpolateFV,interpolateFV_x,interpolateFV_y,interpolateFV_z,gradFV_y, gradFV_z, gradFV_x
        class(sgrid), target, intent(inout) :: this
        real(rkind), dimension(this%nxp, this%nyp, this%nzp, ncnsrv), intent(inout) :: rhs
        real(rkind), dimension(this%nxp, this%nyp, this%nzp), intent(in) :: tauxz,tauyz,tauzz
        real(rkind), dimension(this%nxp, this%nyp, this%nzp), intent(in) :: qz

        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: flux, ddz_p,p_int, w_int, ddz_wp
        real(rkind), dimension(:,:,:), pointer :: ztmp1,ztmp2

        ztmp1 => this%zbuf(:,:,:,1); ztmp2 => this%zbuf(:,:,:,2)

        flux = this%Wcnsrv(:,:,:,mom_index+2)*this%u          - tauxz ! x-momentum
        if(this%use_CnsrvSurfaceTension) then

            flux = flux - this%mix%surfaceTension_fxz

        endif
        call transpose_y_to_z(flux,ztmp1,this%decomp)
        call this%der%ddz(ztmp1,ztmp2,-this%z_bc(1),-this%z_bc(2)) ! Anti-symmetric for all but z-momentum
        call transpose_z_to_y(ztmp2,flux,this%decomp)
        rhs(:,:,:,mom_index  ) = rhs(:,:,:,mom_index  ) - flux
        this%zflux_x = flux

        flux = this%Wcnsrv(:,:,:,mom_index+2)*this%v          - tauyz ! y-momentum
        if(this%use_CnsrvSurfaceTension) then

            flux = flux - this%mix%surfaceTension_fyz

        endif
        call transpose_y_to_z(flux,ztmp1,this%decomp)
        call this%der%ddz(ztmp1,ztmp2,-this%z_bc(1),-this%z_bc(2)) ! Anti-symmetric for all but z-momentum
        call transpose_z_to_y(ztmp2,flux,this%decomp)
        rhs(:,:,:,mom_index+1) = rhs(:,:,:,mom_index+1) - flux
        this%zflux_y = flux

       ! call interpolateFV_z(this%decomp,this%interpMid,this%p,p_int,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
       ! call gradFV_z(this%decomp,this%derStagg,p_int,ddz_p,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
       ! flux = this%Wcnsrv(:,:,:,mom_index+2)*this%w  - tauzz 

        flux = this%Wcnsrv(:,:,:,mom_index+2)*this%w + this%p - tauzz ! z-momentum
        if(this%use_CnsrvSurfaceTension) then

            flux = flux - this%mix%surfaceTension_fzz

        endif
        call transpose_y_to_z(flux,ztmp1,this%decomp)
        call this%der%ddz(ztmp1,ztmp2, this%z_bc(1), this%z_bc(2)) ! Symmetric for z-momentum
        call transpose_z_to_y(ztmp2,flux,this%decomp)
        rhs(:,:,:,mom_index+2) = rhs(:,:,:,mom_index+2) - flux
        this%zflux_z = flux

        !flux = (this%Wcnsrv(:,:,:, TE_index  )  - tauzz)*this%w -this%u*tauxz - this%v*tauyz + qz
        flux = (this%Wcnsrv(:,:,:, TE_index  ) + this%p - tauzz)*this%w - this%u*tauxz - this%v*tauyz + qz ! Total Energy

        if(this%use_CnsrvSurfaceTension) then

            flux = flux - this%mix%surfaceTension_fxz*this%u - this%mix%surfaceTension_fyz*this%v - this%mix%surfaceTension_fzz*this%w

        endif

       ! call interpolateFV_z(this%decomp,this%interpMid,this%w,w_int,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
       ! call gradFV_z(this%decomp,this%derStagg,p_int*w_int,ddz_wp,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        call transpose_y_to_z(flux,ztmp1,this%decomp)
        call this%der%ddz(ztmp1,ztmp2,-this%z_bc(1),-this%z_bc(2)) ! Anti-symmetric for all but z-momentum
        call transpose_z_to_y(ztmp2,flux,this%decomp)
        rhs(:,:,:, TE_index  ) = rhs(:,:,:, TE_index  ) - flux 
        this%zflux_e = flux

    end subroutine

    subroutine filter(this,arr,myfil,numtimes,x_bc_,y_bc_,z_bc_)
        class(sgrid), target, intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(inout) :: arr
        type(filters), target, optional, intent(in) :: myfil
        integer, optional, intent(in) :: numtimes
        integer, dimension(2), optional, intent(in) :: x_bc_, y_bc_, z_bc_
        integer, dimension(2) :: x_bc, y_bc, z_bc
        
        type(filters), pointer :: fil2use
        integer :: times2fil
        real(rkind), dimension(:,:,:), pointer :: tmp_in_y, tmp1_in_x, tmp1_in_z, tmp2_in_x, tmp2_in_z
        integer :: lastx, lasty, lastz, idx


        if (present(myfil)) then
            fil2use => myfil
        else
            fil2use => this%fil
        end if 

        if (present(numtimes)) then
            times2fil = numtimes
        else
            times2fil = 1
        end if

        ! Allocate pointers for the needed buffers 
        ! Atleast 2 buffers in x and z are assumed
        ! Last two buffers are occupied

        lastx = size(this%xbuf,4)
        lasty = size(this%ybuf,4)
        lastz = size(this%zbuf,4)

        tmp1_in_x => this%xbuf(:,:,:,lastx)
        tmp2_in_x => this%xbuf(:,:,:,lastx-1)
        tmp_in_y => this%ybuf(:,:,:,lasty)
        tmp1_in_z => this%zbuf(:,:,:,lastz)
        tmp2_in_z => this%zbuf(:,:,:,lastz-1)

        x_bc = 0; if (present(x_bc_)) x_bc = x_bc_
        y_bc = 0; if (present(y_bc_)) y_bc = y_bc_
        z_bc = 0; if (present(z_bc_)) z_bc = z_bc_
        
        ! First filter in y
        call fil2use%filtery(arr,tmp_in_y,y_bc(1),y_bc(2))
        ! Subsequent refilters 
        do idx = 1,times2fil-1
            arr = tmp_in_y
            call fil2use%filtery(arr,tmp_in_y,y_bc(1),y_bc(2))
        end do
        
        ! Then transpose to x
        call transpose_y_to_x(tmp_in_y,tmp1_in_x,this%decomp)

        ! First filter in x
        call fil2use%filterx(tmp1_in_x,tmp2_in_x,x_bc(1),x_bc(2))
        ! Subsequent refilters
        do idx = 1,times2fil-1
            tmp1_in_x = tmp2_in_x
            call fil2use%filterx(tmp1_in_x,tmp2_in_x,x_bc(1),x_bc(2))
        end do 

        ! Now transpose back to y
        call transpose_x_to_y(tmp2_in_x,tmp_in_y,this%decomp)

        ! Now transpose to z
        call transpose_y_to_z(tmp_in_y,tmp1_in_z,this%decomp)

        !First filter in z
        call fil2use%filterz(tmp1_in_z,tmp2_in_z,z_bc(1),z_bc(2))
        ! Subsequent refilters
        do idx = 1,times2fil-1
            tmp1_in_z = tmp2_in_z
            call fil2use%filterz(tmp1_in_z,tmp2_in_z,z_bc(1),z_bc(2))
        end do 

        ! Now transpose back to y
        call transpose_z_to_y(tmp2_in_z,arr,this%decomp)

        ! Finished
    end subroutine
   
    subroutine getPhysicalProperties(this)
        use exits,      only: GracefulExit
        class(sgrid), intent(inout) :: this

        if (this%mix%ns > 2) then
            call GracefulExit("Number of species must be 1 or 2. for current &
                               implementation of getPhysicalProperties",928)
        endif

        ! If inviscid set everything to zero (otherwise use a model)
        this%mu   = this%phys_mu1   * this%mix%material(1)%VF
        this%bulk = this%phys_bulk1 * this%mix%material(1)%VF
        this%mix%material(1)%physmu = this%phys_mu1
        this%mix%material(2)%physmu = this%phys_mu2

        if (this%mix%ns .eq. 2) then
            this%mu   = this%mu   + this%phys_mu2   * this%mix%material(2)%VF
            this%bulk = this%bulk + this%phys_bulk2 * this%mix%material(2)%VF
        endif

        if (this%PTeqb) then
            this%kap  = this%phys_kap1  * this%mix%material(1)%VF
            if (this%mix%ns .eq. 2) then
                this%kap  = this%kap  + this%phys_kap2  * this%mix%material(2)%VF
            endif
        endif

    end subroutine  

    ! Get tau_ij for momentum equation and simulataneously calculate the viscous work
    ! for the material hydrodynamic equations
    subroutine get_tau(this,duidxj,viscwork)
        class(sgrid), target, intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,9), target, intent(inout) :: duidxj
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),           intent(out)   :: viscwork
        real(rkind), dimension(this%nxp, this%nyp, this%nzp,9), target :: NCbuff, Mubuff   
        real(rkind), dimension(:,:,:), pointer :: d2udx2,d2udy2,d2udz2,d2vdx2,d2vdy2,d2vdz2,d2wdx2,d2wdy2,d2wdz2
        real(rkind), dimension(:,:,:), pointer :: dmudx,dmudy,dmudz
        real(rkind), dimension(:,:,:), pointer :: dbulkdx,dbulkdy,dbulkdz
        real(rkind), dimension(:,:,:), pointer :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
        real(rkind), dimension(:,:,:), pointer :: lambda, bambda

        lambda => this%ybuf(:,:,:,1)
        bambda => this%ybuf(:,:,:,2)

        dudx => duidxj(:,:,:,1); dudy => duidxj(:,:,:,2); dudz => duidxj(:,:,:,3);
        dvdx => duidxj(:,:,:,4); dvdy => duidxj(:,:,:,5); dvdz => duidxj(:,:,:,6);
        dwdx => duidxj(:,:,:,7); dwdy => duidxj(:,:,:,8); dwdz => duidxj(:,:,:,9);
       
        ! Compute the multiplying factors (thermo-shit)
        bambda = (four/three)*this%mu + this%bulk
        lambda = this%bulk - (two/three)*this%mu
    
    
        ! Step 1: Get tau_12  (dudy is destroyed)
        dudy =  dudy + dvdx
        viscwork  = (this%mu*dudy) * (dudy)  ! tau_12 * (S_12+S_21) (Since symmetric)
        dudy = this%mu*dudy
        !tauxyidz = 2
    
        ! Step 2: Get tau_13 (dudz is destroyed)
        dudz = dudz + dwdx
        viscwork  = viscwork + (this%mu*dudz) * (dudz)  ! tau_13 *(S_13 + S_31)
        dudz = this%mu*dudz
        !tauxzidx = 3

        ! Step 3: Get tau_23 (dvdz is destroyed)
        dvdz = dvdz + dwdy
        viscwork  = viscwork + (this%mu*dvdz) * (half*dvdz)  ! tau_23 * (S_23 + S_32)
        dvdz = this%mu*dvdz
        !tauyzidx = 6

        ! Step 4: Get tau_11 (dvdx is destroyed)
        dvdx = bambda*dudx + lambda*(dvdy + dwdz)
        viscwork  = viscwork + dvdx * dudx  ! tau_11 * S_11
        !tauxxidx = 4

        ! Step 5: Get tau_22 (dwdx is destroyed)
        dwdx = bambda*dvdy + lambda*(dudx + dwdz)
        viscwork  = viscwork + dwdx * dvdy  ! tau_22 * S_22
        !tauyyidx = 7

        ! Step 6: Get tau_33 (dwdy is destroyed)
        dwdy = bambda*dwdz + lambda*(dudx + dvdy)
        viscwork  = viscwork + dwdy * dwdz  ! tau_33 * S_33
        !tauzzidx = 8
    
        ! Done 
    end subroutine 


    subroutine get_tauSum(this,tauSum)
        use operators, only: divergence,gradient,divergenceFV,interpolateFV,interpolateFV_x,interpolateFV_y,interpolateFV_z,gradFV_x, gradFV_y, gradFV_z
        class(sgrid), target, intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(inout) :: tauSum
        real(rkind), dimension(this%nxp, this%nyp, this%nzp) :: dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz
        real(rkind), dimension(this%nxp, this%nyp, this%nzp,3) :: w_int, v_int, u_int
        real(rkind), dimension(:,:,:), pointer :: lambda, bambda

        lambda => this%ybuf(:,:,:,1)
        bambda => this%ybuf(:,:,:,2)

        !call gradient(this%decomp,this%derCD06,this%u, dudx, dudy, dudz,-this%x_bc,  this%y_bc,this%z_bc)
        !call gradient(this%decomp,this%derCD06,this%v, dvdx, dvdy, dvdz,this%x_bc, -this%y_bc,this%z_bc)
        !call gradient(this%decomp,this%derCD06,this%w, dwdx, dwdy, dwdz,this%x_bc,  this%y_bc,-this%z_bc)

        call interpolateFV_x(this%decomp,this%interpMid,this%w,w_int(:,:,:,1),this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        call interpolateFV_y(this%decomp,this%interpMid,this%w,w_int(:,:,:,2),this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        call interpolateFV_z(this%decomp,this%interpMid,this%w,w_int(:,:,:,3),this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)

        call interpolateFV_x(this%decomp,this%interpMid,this%v,v_int(:,:,:,1),this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        call interpolateFV_y(this%decomp,this%interpMid,this%v,v_int(:,:,:,2),this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        call interpolateFV_z(this%decomp,this%interpMid,this%v,v_int(:,:,:,3),this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)

        call interpolateFV_x(this%decomp,this%interpMid,this%u,u_int(:,:,:,1),this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        call interpolateFV_y(this%decomp,this%interpMid,this%u,u_int(:,:,:,2),this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        call interpolateFV_z(this%decomp,this%interpMid,this%u,u_int(:,:,:,3),this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)

        call gradFV_x(this%decomp,this%derStagg,u_int(:,:,:,1),dudx,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        call gradFV_y(this%decomp,this%derStagg,u_int(:,:,:,2),dudy,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        call gradFV_z(this%decomp,this%derStagg,u_int(:,:,:,3),dudz,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)

        call gradFV_x(this%decomp,this%derStagg,v_int(:,:,:,1),dvdx,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        call gradFV_y(this%decomp,this%derStagg,v_int(:,:,:,2),dvdy,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        call gradFV_z(this%decomp,this%derStagg,v_int(:,:,:,3),dvdz,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)

        call gradFV_x(this%decomp,this%derStagg,w_int(:,:,:,1),dwdx,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        call gradFV_y(this%decomp,this%derStagg,w_int(:,:,:,2),dwdy,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        call gradFV_z(this%decomp,this%derStagg,w_int(:,:,:,3),dwdz,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)

        this%dudx = dudx; this%dudy = dudy; this%dudz = dudz;
        this%dwdx = dwdx; this%dwdy = dwdy; this%dwdz = dwdz;
        this%dvdx = dvdx; this%dvdy = dvdy; this%dvdz = dvdz;

        ! Compute the multiplying factors (thermo-shit)
        bambda = (four/three)*this%mu + this%bulk
        lambda = this%bulk - (two/three)*this%mu


        ! Step 1: Get tau_12  (dudy is destroyed)
        this%tauxy =  this%mu*(dudy + dvdx)
        !tauxyidz = 2

        ! Step 2: Get tau_13 (dudz is destroyed)
        this%tauxz = this%mu*(dudz + dwdx)
        !tauxzidx = 3

        ! Step 3: Get tau_23 (dvdz is destroyed)
        this%tauyz = this%mu*(dvdz + dwdy)
        !tauyzidx = 6

        ! Step 4: Get tau_11 (dvdx is destroyed)
        this%tauxx = bambda*dudx + lambda*(dvdy + dwdz)

        ! Step 5: Get tau_22 (dwdx is destroyed)
        this%tauyy = bambda*dvdy + lambda*(dudx + dwdz)
        !tauyyidx = 7

        ! Step 6: Get tau_33 (dwdy is destroyed)
        this%tauzz = bambda*dwdz + lambda*(dudx + dvdy)
        !tauzzidx = 8

        tauSum =  this%tauxx*dudx + this%tauyy*dvdy + this%tauzz*dwdz &
                + this%tauxy*dudy + this%tauxz*dudz + this%tauyz*dvdz &
                + this%tauxy*dvdx + this%tauxz*dwdx + this%tauyz*dwdy

        this%tauSum = tauSum
        ! Done 
    end subroutine

    subroutine get_tauStagg(this,duidxj,duidxj_int,duidxj_s)
        use operators, only : interpolateFV_x, interpolateFV_y, interpolateFV_z
        class(sgrid), target, intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,9), target, intent(inout) :: duidxj, duidxj_s
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,12), target,intent(inout) :: duidxj_int
        real(rkind), dimension(this%nxp, this%nyp, this%nzp,9), target :: NCbuff, Mubuff
        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: mu_x, mu_y, mu_z,bulk_x, bulk_y, bulk_z
        real(rkind), dimension(:,:,:), pointer :: dmudx,dmudy,dmudz
        real(rkind), dimension(:,:,:), pointer :: dbulkdx,dbulkdy,dbulkdz
        real(rkind), dimension(:,:,:), pointer :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
        real(rkind), dimension(:,:,:), pointer :: lambda, bambda
        real(rkind), dimension(:,:,:), pointer :: dudx_s,dudy_s,dudz_s,dvdx_s,dvdy_s,dvdz_s,dwdx_s,dwdy_s,dwdz_s
        real(rkind), dimension(:,:,:), pointer :: dvdy_x,dwdz_x, dvdx_y, dwdx_z
        real(rkind), dimension(:,:,:), pointer :: dudx_y, dwdz_y, dudy_x, dwdy_z
        real(rkind), dimension(:,:,:), pointer :: dudx_z, dvdy_z, dudz_x, dvdz_y

        call interpolateFV_x(this%decomp,this%interpMid,this%mu,mu_x,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        call interpolateFV_y(this%decomp,this%interpMid,this%mu,mu_y,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        call interpolateFV_z(this%decomp,this%interpMid,this%mu,mu_z,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)

        call interpolateFV_x(this%decomp,this%interpMid,this%bulk,bulk_x,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        call interpolateFV_y(this%decomp,this%interpMid,this%bulk,bulk_y,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        call interpolateFV_z(this%decomp,this%interpMid,this%bulk,bulk_z,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)

        lambda => this%ybuf(:,:,:,1)
        bambda => this%ybuf(:,:,:,2)

        dudx => duidxj(:,:,:,1); dudy => duidxj(:,:,:,2); dudz => duidxj(:,:,:,3);
        dvdx => duidxj(:,:,:,4); dvdy => duidxj(:,:,:,5); dvdz => duidxj(:,:,:,6);
        dwdx => duidxj(:,:,:,7); dwdy => duidxj(:,:,:,8); dwdz => duidxj(:,:,:,9);

        dvdy_x => duidxj_int(:,:,:,1); dudx_y => duidxj_int(:,:,:,2); dudx_z => duidxj_int(:,:,:,3);
        dwdz_x => duidxj_int(:,:,:,4); dwdz_y => duidxj_int(:,:,:,5); dvdy_z => duidxj_int(:,:,:,6);
        dvdx_y => duidxj_int(:,:,:,7); dudy_x => duidxj_int(:,:,:,8); dudz_x => duidxj_int(:,:,:,9);
        dwdx_z => duidxj_int(:,:,:,10); dwdy_z => duidxj_int(:,:,:,11); dvdz_y => duidxj_int(:,:,:,12);

        dudx_s => duidxj_s(:,:,:,1); dudy_s => duidxj_s(:,:,:,2); dudz_s => duidxj_s(:,:,:,3);
        dvdx_s => duidxj_s(:,:,:,4); dvdy_s => duidxj_s(:,:,:,5); dvdz_s => duidxj_s(:,:,:,6);
        dwdx_s => duidxj_s(:,:,:,7); dwdy_s => duidxj_s(:,:,:,8); dwdz_s => duidxj_s(:,:,:,9);

        ! Compute the multiplying factors (thermo-shit)
        bambda = (four/three)*mu_x + bulk_x
        lambda = bulk_x - (two/three)*mu_x


        ! Step 1: Get tau_12  (dudy is destroyed)
        dudy =  dudy_x + dvdx_s
        dudy = mu_x*dudy
        !tauxyidz = 2

        ! Step 2: Get tau_13 (dudz is destroyed)
        dudz = dudz_x + dwdx_s
        dudz = mu_x*dudz
        !tauxzidx = 3

        ! Step 3: Get tau_23 (dvdz is destroyed)
        dvdz = dvdz_y +  dwdy_s
        dvdz = mu_y*dvdz
        !tauyzidx = 6

        ! Step 4: Get tau_11 (dvdx is destroyed)
        dvdx = bambda*dudx_s + lambda*(dvdy_x + dwdz_x)
        !tauxxidx = 4

        bambda = (four/three)*mu_y + bulk_y
        lambda = bulk_y - (two/three)*mu_y

        ! Step 5: Get tau_22 (dwdx is destroyed)
        dwdx = bambda*dvdy_s + lambda*(dudx_y + dwdz_y)
        !tauyyidx = 7
     
        bambda = (four/three)*mu_z + bulk_z
        lambda = bulk_z - (two/three)*mu_z

        ! Step 6: Get tau_33 (dwdy is destroyed)
        dwdy = bambda*dwdz_s + lambda*(dudx_z + dvdy_z)
        !tauzzidx = 8

        ! tau 21
        dvdy_x = mu_y*(dvdx_y + dudy_s)
        !tauyxid = 1
    
        ! tau 31
        dudx_y  = mu_z*( dwdx_z + dudz_s )

        ! tau 32 
        dudx_z = mu_z*(dvdz_s + dwdy_z)
        ! Done 
    end subroutine

    subroutine get_q(this,qx,qy,qz)
        use exits, only: nancheck
        class(sgrid), target, intent(inout) :: this
        real(rkind), dimension(this%nxp, this%nyp, this%nzp), intent(inout) :: qx,qy,qz

        ! integer :: i
        real(rkind), dimension(:,:,:), pointer :: tmp1_in_x, tmp2_in_x, tmp1_in_y, tmp1_in_z, tmp2_in_z
        type(derivatives), pointer :: der

        der => this%der

        tmp1_in_x => this%xbuf(:,:,:,1)
        tmp2_in_x => this%xbuf(:,:,:,2)

        tmp1_in_z => this%zbuf(:,:,:,1)
        tmp2_in_z => this%zbuf(:,:,:,2)

        tmp1_in_y => this%ybuf(:,:,:,1)

        ! Species enthalpy diffusion is computed earlier in SolidMixture

        ! Step 1: Get qy (dvdy is destroyed)
        call der%ddy(this%T,tmp1_in_y,this%y_bc(1),this%y_bc(2))
        qy = qy - this%kap*tmp1_in_y

        ! Step 2: Get qx (dudx is destroyed)
        call transpose_y_to_x(this%T,tmp1_in_x,this%decomp)
        call der%ddx(tmp1_in_x,tmp2_in_x,this%x_bc(1),this%x_bc(2))
        call transpose_x_to_y(tmp2_in_x,tmp1_in_y,this%decomp)
        qx = qx - this%kap*tmp1_in_y

        ! Step 3: Get qz (dwdz is destroyed)
        call transpose_y_to_z(this%T,tmp1_in_z,this%decomp)
        call der%ddz(tmp1_in_z,tmp2_in_z,this%z_bc(1),this%z_bc(2))
        call transpose_z_to_y(tmp2_in_z,tmp1_in_y)
        qz = qz - this%kap*tmp1_in_y

        ! Done
    end subroutine

    subroutine get_qLAD(this,qx,qy,qz,qDiv)
        use exits, only: nancheck
        use operators, only: divergenceFV, interpolateFV,interpolateFV_x, interpolateFV_y, interpolateFV_z, gradFV_N2Fx, gradFV_N2Fy,gradFV_N2Fz

        class(sgrid), target, intent(inout) :: this
        real(rkind), dimension(this%nxp, this%nyp, this%nzp), intent(inout) :: qx,qy,qz,qDiv

        ! integer :: i
        real(rkind), dimension(:,:,:), pointer :: tmp1_in_x, tmp2_in_x,tmp1_in_y, tmp1_in_z, tmp2_in_z
        type(derivatives), pointer :: der
        real(rkind), dimension(this%nxp, this%nyp, this%nzp,3) :: kap_int

        qx = 0
        qy = 0
        qz = 0

        call interpolateFV_x(this%decomp,this%interpMid,this%kap,kap_int(:,:,:,1),this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        call interpolateFV_y(this%decomp,this%interpMid,this%kap,kap_int(:,:,:,2),this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        call interpolateFV_z(this%decomp,this%interpMid,this%kap,kap_int(:,:,:,3),this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)

        call gradFV_N2Fx(this%decomp,this%derStagg,this%T,qx,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        call gradFV_N2Fy(this%decomp,this%derStagg,this%T,qy,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        call gradFV_N2Fz(this%decomp,this%derStagg,this%T,qz,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)

        qx = -kap_int(:,:,:,1)*qx
        qy = -kap_int(:,:,:,2)*qy
        qz = -kap_int(:,:,:,3)*qz

        call divergenceFV(this%decomp,this%derStagg,qx,qy,qz,qDiv,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        ! Done
    end subroutine

end module 
