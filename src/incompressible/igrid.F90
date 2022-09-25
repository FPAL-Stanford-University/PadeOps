module IncompressibleGrid
    use kind_parameters, only: rkind, clen
    use constants, only: imi, zero,one,two,three,half,fourth, pi, kappa 
    use GridMod, only: grid
    use gridtools, only: alloc_buffs, destroy_buffs
    use igrid_hooks!, only: setDirichletBC_Temp, set_Reference_Temperature, meshgen_WallM, initfields_wallM, set_planes_io, set_KS_planes_io 
    use decomp_2d
    use StaggOpsMod, only: staggOps  
    use exits, only: GracefulExit, message, check_exit
    use spectralMod, only: spectral  
    !use PoissonMod, only: poisson
    use mpi 
    use reductions, only: p_maxval, p_sum, p_minval
    use timer, only: tic, toc
    use PadePoissonMod, only: Padepoisson 
    use sgsmod_igrid, only: sgs_igrid
    use numerics
    !use cd06staggstuff, only: cd06stagg
    use cf90stuff, only: cf90
    use TurbineMod, only: TurbineArray 
    use kspreprocessing, only: ksprep  
    use PadeDerOps, only: Pade6Stagg
    use Fringemethod, only: fringe
    use angleControl, only: angCont
    use forcingmod,   only: HIT_shell_forcing
    use scalar_igridMod, only: scalar_igrid 
    use io_hdf5_stuff, only: io_hdf5 
    use PoissonPeriodicMod, only: PoissonPeriodic

    implicit none

    external :: MPI_BCAST, MPI_RECV, MPI_SEND, MPI_REDUCE

    private
    public :: igrid, wBC_bottom, wBC_top  

    complex(rkind), parameter :: zeroC = zero + imi*zero 

    ! Allow non-zero value (isEven) 
    integer :: ierr 
    integer :: topWall = 1, botWall = 1, botBC_Temp = 1, topBC_Temp = 0

    !! BC convention: 

    !! +1: even extension
    !! -1: odd extension
    !!  0: sided stencil
    !!  Default values: 
    integer :: wBC_bottom     = -1, wBC_top     = -1
    integer :: uBC_bottom     =  0, uBC_top     =  1
    integer :: vBC_bottom     =  0, vBC_top     =  1
    integer :: TBC_bottom     =  1, TBC_top     =  0
    integer :: WdUdzBC_bottom = -1, WdUdzBC_top =  1
    integer :: WdVdzBC_bottom = -1, WdVdzBC_top =  1
    integer :: WdWdzBC_bottom = -1, WdWdzBC_top = -1
    integer :: WWBC_bottom    =  1, WWBC_top    =  1
    integer :: UWBC_bottom    = -1, UWBC_top    = -1
    integer :: VWBC_bottom    = -1, VWBC_top    = -1
    integer :: WTBC_bottom    = -1, WTBC_top    = -1
    integer :: dUdzBC_bottom  =  0, dUdzBC_top  =  -1
    integer :: dVdzBC_bottom  =  0, dVdzBC_top  =  -1
    integer :: dWdzBC_bottom  =  1, dWdzBC_top  =  1
    integer :: dTdzBC_bottom  =  -1, dTdzBC_top  =  0
    integer :: WdTdzBC_bottom =   1, WdTdzBC_top = 0

    type :: igrid
        
        character(clen) :: inputDir
        integer :: headerfid = 12345   

        ! Variables common to grid
        integer :: nx, ny, nz, t_datadump, t_restartdump
        real(rkind) :: dt, tstop, CFL, CviscDT, dx, dy, dz, tsim
        character(len=clen) ::  outputdir
        real(rkind), dimension(:,:,:,:), allocatable :: mesh
        real(rkind) :: zBot, zTop, zMid
        character(len=clen) :: filter_x          ! What filter to use in X: "cf90", "gaussian", "lstsq", "spectral"
        character(len=clen) :: filter_y          ! What filter to use in X: "cf90", "gaussian", "lstsq", "spectral" 
        character(len=clen) :: filter_z          ! What filter to use in X: "cf90", "gaussian", "lstsq", "spectral" 
        integer      :: step, nsteps = 999999
        logical :: Am_I_Primary = .true.

        type(decomp_info), allocatable :: gpC, gpE
        type(decomp_info), pointer :: Sp_gpC, Sp_gpE
        type(spectral), allocatable :: spectE, spectC
        type(staggOps), allocatable :: Ops, OpsPP
        type(sgs_igrid), allocatable :: sgsmodel

        real(rkind), dimension(:,:,:,:), allocatable :: PfieldsC
        real(rkind), dimension(:,:,:,:), allocatable :: PfieldsE
        !type(cd06stagg), allocatable :: derW, derWW, derSO, derSE, derT
        type(Pade6Stagg), allocatable :: Pade6opZ
        type(cf90),      allocatable :: filzE, filzC

        complex(rkind), dimension(:,:,:,:), allocatable :: SfieldsC
        complex(rkind), dimension(:,:,:,:), allocatable :: SfieldsE


        type(padepoisson), allocatable :: padepoiss
        real(rkind), dimension(:,:,:), allocatable :: divergence

        real(rkind), dimension(:,:,:), pointer :: u, v, wC, w, uE, vE, T, TE
        complex(rkind), dimension(:,:,:), pointer :: uhat, vhat, whatC, what, That, TEhat, uEhat, vEhat
        
        complex(rkind), dimension(:,:,:), pointer :: uhat1, vhat1, what1, That1
        complex(rkind), dimension(:,:,:), pointer :: uhat2, vhat2, what2, That2
        complex(rkind), dimension(:,:,:), pointer :: uhat3, vhat3, what3, That3
        complex(rkind), dimension(:,:,:), pointer :: uhat4, vhat4, what4, That4
        complex(rkind), dimension(:,:,:,:), allocatable :: SfieldsC2, SfieldsE2
        complex(rkind), dimension(:,:,:,:), allocatable :: uExtra, vExtra, wExtra, TExtra 
        complex(rkind), dimension(:,:,:,:), allocatable :: uRHSExtra, vRHSExtra, wRHSExtra, TRHSExtra 

        real(rkind), dimension(:,:,:), pointer :: ox,oy,oz
        complex(rkind), dimension(:,:,:), pointer :: T_rhs, T_Orhs

        complex(rkind), dimension(:,:,:), allocatable :: uBase, Tbase, vBase, dTdxH, dTdyH, dTdzH, dTdzHC
        real(rkind), dimension(:,:,:), allocatable :: dTdxC, dTdyC, dTdzE, dTdzC, dTdxE, dTdyE


        real(rkind), dimension(:,:,:,:), allocatable, public :: rbuffxC, rbuffyC, rbuffzC
        real(rkind), dimension(:,:,:,:), allocatable :: rbuffxE, rbuffyE, rbuffzE
        
        complex(rkind), dimension(:,:,:,:), allocatable :: cbuffyC, cbuffzC!, cbuffxC
        complex(rkind), dimension(:,:,:,:), allocatable :: cbuffyE, cbuffzE

        complex(rkind), dimension(:,:,:,:), allocatable :: rhsC, rhsE, OrhsC, OrhsE 
        real(rkind), dimension(:,:,:,:), allocatable :: duidxjC, duidxjE 
        complex(rkind), dimension(:,:,:,:), allocatable :: duidxjChat, duidxjEhat
        complex(rkind), dimension(:,:,:), allocatable :: d2udz2hatC, d2vdz2hatC,d2wdz2hatE, d2Tdz2hatC
        complex(rkind), dimension(:,:,:), pointer:: u_rhs, v_rhs, wC_rhs, w_rhs 
        complex(rkind), dimension(:,:,:), pointer:: u_Orhs, v_Orhs, w_Orhs

        real(rkind), dimension(:,:,:), allocatable :: rDampC, rDampE         
        real(rkind) :: Re, G_Geostrophic, G_alpha, frameAngle, dtby2, meanfact, Tref, dPfdx, dPfdy, dPfdz
        complex(rkind), dimension(:,:,:), allocatable :: GxHat, GyHat, GxHat_Edge, GyHat_Edge
        real(rkind) :: Ro = 1.d5, Fr = 1000.d0, PrandtlFluid = 1.d0, BulkRichardson = 0.d90
        logical :: assume_fplane = .true.
        real(rkind) :: coriolis_omegaY, coriolis_omegaZ, coriolis_omegaX 
        integer :: nxZ, nyZ
   
        character(len=clen) :: dtlimit
        integer :: BuoyancyTermType = 0 
        real(rkind) :: Ra, BuoyancyFact = 0.d0

        integer :: moistureIndex = 1
        real(rkind) :: moistureFactor = 0.61d0 ! converts g/kg to K

        logical :: periodicInZ  = .false. 
        logical :: newTimeStep = .true., computeVorticity = .false.  
        integer :: timeSteppingScheme = 0 
        integer :: runID, t_start_planeDump, t_stop_planeDump, t_planeDump, t_DivergenceCheck
        integer :: t_start_pointProbe, t_stop_pointProbe, t_pointProbe
        logical :: useCoriolis = .true. , isStratified = .false., useSponge = .false., useMoisture = .false.
        logical :: useExtraForcing = .false., useGeostrophicForcing = .false., isInviscid = .false.  
        logical :: useSGS = .false., computeTurbinePressure = .false.  
        logical :: UseDealiasFilterVert = .false.
        logical :: useDynamicProcedure 
        logical :: useCFL = .false., donot_dealias = .false.   
        logical :: dumpPlanes = .false., useWindTurbines = .false. 

        complex(rkind), dimension(:,:,:), allocatable :: dPf_dxhat

        real(rkind) :: latitude, max_nuSGS, invObLength, Tsurf0, Tsurf, dTsurf_dt, ThetaRef, TsurfTop, dTsurfTop_dt

        real(rkind) :: dtOld, dtRat, Tmn, wTh_surf
        integer :: wallMType, botBC_Temp, topBC_Temp 

        ! Statistics to compute 
        real(rkind), dimension(:), allocatable :: runningSum_sc, inst_horz_avg, runningSum_sc_turb, runningSum_turb, inst_horz_avg_turb, debugavg, debuginst
        real(rkind), dimension(:,:), allocatable :: zStats2dump, runningSum, TemporalMnNOW, horzavgstats
        real(rkind), dimension(:,:,:,:), allocatable :: stats3D
        real(rkind), dimension(:,:,:),   allocatable :: xspectra_mean
        !---pointers for horizontally- and time-averaged statistics; linked to horzavgstats------
        real(rkind), dimension(:), pointer :: u_mean, v_mean, w_mean, uu_mean, uv_mean, uw_mean, vv_mean, vw_mean, ww_mean, disperuw_mean, dispervw_mean
        real(rkind), dimension(:), pointer :: mkeadv_mean, mkett_mean, mkedisp_mean, tkeadv_mean, tkett_mean, tkeprod_mean
        real(rkind), dimension(:), pointer :: TT_mean, wT_mean, vT_mean, uT_mean, T_mean
        real(rkind), dimension(:), pointer :: p_mean, mkept_mean, tkept_mean
        real(rkind), dimension(:), pointer :: mkevdif_mean, mkevdsp_mean, tkevdif_mean, tkevdsp_mean
        real(rkind), dimension(:), pointer :: mkesgst_mean, mkesgsd_mean, tkesgst_mean, tkesgsd_mean
        real(rkind), dimension(:), pointer :: tau11_mean, tau12_mean, tau13_mean, tau22_mean, tau23_mean, tau33_mean
        real(rkind), dimension(:), pointer :: q1_mean, q2_mean, q3_mean
        real(rkind), dimension(:), pointer :: turbfx_mean, turbfy_mean, turbfz_mean, tketurbf_mean, mketurbf_mean
        !---pointers for horizontally- and time-averaged statistics; linked to horzavgstats------
        !---pointers for time-averaged statistics; linked to stats3D------
        real(rkind), dimension(:,:,:), pointer :: u_mean3D, v_mean3D, w_mean3D, uu_mean3D, uv_mean3D, uw_mean3D, vv_mean3D, vw_mean3D, ww_mean3D, tketurbtranspx_mean3D, tketurbtranspy_mean3D, tketurbtranspz_mean3D
        real(rkind), dimension(:,:,:), pointer :: TT_mean3D, wT_mean3D, vT_mean3D, uT_mean3D, T_mean3D
        real(rkind), dimension(:,:,:), pointer :: p_mean3D, pw_mean3D, pv_mean3D, pu_mean3D
        real(rkind), dimension(:,:,:), pointer :: viscdisp_mean3D, Siju1_mean3D, Siju2_mean3D, Siju3_mean3D
        real(rkind), dimension(:,:,:), pointer :: tau11_mean3D, tau12_mean3D, tau13_mean3D, tau22_mean3D, tau23_mean3D, tau33_mean3D
        real(rkind), dimension(:,:,:), pointer :: sgsdissp_mean3D, tauu1_mean3D, tauu2_mean3D, tauu3_mean3D, q1_mean3D, q2_mean3D, q3_mean3D
        real(rkind), dimension(:,:,:), pointer :: turbfx_mean3D, turbfy_mean3D, turbfz_mean3D, uturbf_mean3D
        !---pointers for time-averaged statistics; linked to stats3D------
        integer :: tidSUM, tid_StatsDump, tid_compStats, tprev2, tprev1
        logical :: normByustar
        real(rkind) :: tSimStartStats
        real(rkind), dimension(:,:,:,:), allocatable :: F_rhs_sgs
        logical :: computeForcingTerm = .false., deleteInstructions 

        ! Pointers linked to SGS stuff
        real(rkind), dimension(:,:,:,:), pointer :: tauSGS_ij
        real(rkind), dimension(:,:,:)  , pointer :: kappaSGS, nu_SGS, tau13, tau23, kappa_bounding
        real(rkind), dimension(:,:,:)  , pointer :: c_SGS, q1, q2, q3 
        real(rkind), dimension(:,:,:), allocatable :: fbody_x, fbody_y, fbody_z
        logical                                      :: storeFbody

        ! Wind Turbine stuff 
        type(turbineArray), allocatable :: WindTurbineArr

        ! KS preprocessor 
        type(ksprep), allocatable :: LES2KS
        character(len=clen) :: KSoutputdir
        logical :: PreProcessForKS, KSupdated = .false. 
        integer, dimension(:), allocatable :: planes2dumpC_KS, planes2dumpF_KS
        integer :: t_dumpKSprep, KSinitType
        real(rkind), dimension(:,:,:), pointer :: uFil4KS, vFil4KS, wFil4KS
        real(rkind) :: turbPr, KSFilFact 
        real(rkind), dimension(:,:,:), allocatable :: KS_probe_data


        ! Pressure Solver
        logical :: StorePressure = .false., fastCalcPressure = .true. 
        integer :: P_dumpFreq = 10, P_compFreq = 10
        logical :: AlreadyHaveRHS = .false.
        logical :: ComputeDNSPressure = .false., ComputeFringePressure = .false.  
        real(rkind), dimension(:,:,:), allocatable :: pressure, pressure_dns, pressure_fringe, pressure_turbine
        complex(rkind), dimension(:,:,:), allocatable :: urhs_dns,vrhs_dns,wrhs_dns,urhs_fringe,vrhs_fringe,wrhs_fringe
        complex(rkind), dimension(:,:,:), allocatable :: urhs_turbine, vrhs_turbine, wrhs_turbine

        logical :: Dump_NU_SGS = .false., Dump_KAPPA_SGS = .false. 

        ! Rapid and slow decomposition
        real(rkind), dimension(:,:,:), allocatable :: prapid, uM, vM, wM, pslow
        real(rkind), dimension(:,:,:), allocatable :: dumdx, dumdy, dumdz
        real(rkind), dimension(:,:,:), allocatable :: dvmdx, dvmdy, dvmdz
        real(rkind), dimension(:,:,:), allocatable :: dwmdx, dwmdy, dwmdz
        logical :: computeRapidSlowPressure 
        type(PoissonPeriodic) :: poiss_periodic

        ! Stats
        logical :: timeAvgFullFields, computeSpectra

        ! System Interactions 
        logical :: useSystemInteractions = .false. 
        integer :: tSystemInteractions = 10
        character(len=clen) :: controlDir = "null"

        integer, dimension(:), allocatable :: xplanes, yplanes, zplanes
        ! Note that c_SGS is linked to a variable that is constant along & 
        ! i, j but is still stored as a full 3 rank array. This is mostly done to
        ! make it convenient us to later do transposes or to compute Sij.

        ! Probes
        logical :: useProbes = .false., doIhaveAnyProbes = .false.  
        integer, dimension(:,:), allocatable :: probes
        integer :: nprobes, probeTimeLimit = 1000000, probeStartStep = 0 
        real(rkind), dimension(:,:,:), allocatable :: probe_data
        integer :: tpro

        ! Spin up for initialization
        logical :: InitSpinUp = .false. 
        real(rkind) :: Tstop_InitSpinUp

        ! Fringe
        logical                           :: useFringe = .false., usedoublefringex = .false. 
        type(fringe), allocatable, public :: fringe_x, fringe_x1, fringe_x2

        ! Control
        logical                            :: useControl = .false.
        type(angCont), allocatable, public :: angCont_yaw
        real(rkind) :: angleHubHeight, totalAngle, wFilt, restartPhi, deltaGalpha, angleTrigger
        integer :: zHubIndex = 16

        ! HIT Forcing
        logical :: useHITForcing = .false., useforcedStratification = .false.
        logical :: useHITRealSpaceLinearForcing = .false.
        type(HIT_shell_forcing), allocatable :: hitforce
        real(rkind) :: HITForceTimeScale 

        ! Scalars
        logical :: useScalars = .false. 
        type(scalar_igrid), dimension(:), allocatable :: scalars
        integer :: n_scalars = 1

        ! Dump schedule
        integer :: vizDump_Schedule = 0
        real(rkind) :: deltaT_dump, t_NextDump
        logical :: DumpThisStep = .false.

        ! HDF5 IO
        integer :: ioType = 0
        type(io_hdf5) :: viz_hdf5

        logical :: WriteTurbineForce = .false. 
        
        ! budgets on the fly
        logical :: StoreForBudgets = .false. 
        complex(rkind), dimension(:,:,:), pointer :: ucon, vcon, wcon, usgs, vsgs, wsgs, uvisc, vvisc, wvisc, px, py, pz, wb, ucor, vcor, wcor, uturb, pxdns, pydns, pzdns, vturb, wturb, HITforcing_x, HITforcing_y, HITforcing_z 

        integer :: buoyancyDirection 

        contains
            procedure          :: init
            procedure          :: destroy
            procedure          :: printDivergence 
            procedure          :: getMaxKE
            procedure          :: getMeanKE
            procedure          :: getMeanuu
            procedure          :: getMeanuv
            procedure          :: getMeanuw
            procedure          :: getMeanvv
            procedure          :: getMeanvw
            procedure          :: getMeanww
            procedure          :: timeAdvance
            procedure          :: start_io
            procedure          :: finalize_io
            procedure          :: get_dt
            procedure          :: interpolate_cellField_to_edgeField 
            !procedure, private :: init_stats
            procedure, private :: init_stats3D
            procedure, private :: AdamsBashforth
            procedure, private :: TVD_RK3
            procedure, private :: SSP_RK45
            procedure, private :: ComputePressure
            procedure, private :: interp_primitiveVars
            procedure, private :: compute_duidxj
            procedure, private :: compute_Sijmean
            procedure, private :: compute_dTdxi
            procedure, private :: addNonLinearTerm_Rot
            procedure, private :: addBuoyancyTerm
            procedure, private :: addViscousTerm
            procedure, private :: addCoriolisTerm
            procedure, private :: addSponge
            procedure, private :: addExtraForcingTerm 
            procedure, private :: addForcedStratification
            procedure, private :: dumpRestartFile
            procedure, private :: readRestartFile
            procedure, private :: compute_z_mean 
            procedure, private :: compute_deltaT
            procedure, private :: dump_pointProbes
            procedure, private :: dump_stats3D
            procedure, private :: compute_stats3D 
            procedure, private :: dealiasFields
            procedure, private :: dealias_rhs
            procedure, private :: ApplyCompactFilter 
            procedure, private :: addNonLinearTerm_skewSymm
            procedure, private :: populate_rhs
            procedure, private :: populate_rhs_for_budgets
            procedure, private :: populate_RHS_extraTerms 
            procedure, private :: project_and_prep
            procedure          :: wrapup_timestep
            procedure, private :: reset_pointers
            procedure, private :: compute_vorticity
            procedure, private :: finalize_stats3D
            procedure, private :: dump_planes
            procedure, private :: dealiasRealField_C
            procedure          :: dumpFullField 
            procedure          :: dumpSpectralField 
            procedure, private :: dump_scalar_fields
            procedure, private :: dumpVisualizationInfo
            procedure, private :: DeletePrevStats3DFiles
            procedure, private :: Delete_file_if_present
            procedure, private :: updateProbes 
            procedure, private :: dumpProbes
            procedure, private :: correctPressureRotationalForm
            procedure, private :: initialize_scalar_for_InitSpinUp
            procedure, private :: initialize_Rapid_Slow_Pressure_Split
            procedure, private :: compute_RapidSlowPressure_Split
            procedure, private :: dump_visualization_files
            procedure, private :: append_visualization_info
            procedure, private :: initialize_hdf5_io
            procedure, private :: destroy_hdf5_io
            procedure          :: get_geostrophic_forcing
            procedure          :: instrumentForBudgets
            procedure          :: instrumentForBudgets_timeAvg
            procedure          :: instrumentForBudgets_volAvg
            procedure          :: getMomentumTerms
            procedure          :: set_budget_rhs_to_zero
            procedure, private :: advance_SSP_RK45_all_stages
            procedure          :: advance_SSP_RK45_Stage_1 
            procedure          :: advance_SSP_RK45_Stage_2 
            procedure          :: advance_SSP_RK45_Stage_3 
            procedure          :: advance_SSP_RK45_Stage_4 
            procedure          :: advance_SSP_RK45_Stage_5
            procedure          :: getMaxOddballModes 
   end type

contains 

#include "igrid_files/io_stuff.F90"
#include "igrid_files/stats_stuff.F90"
#include "igrid_files/rhs_stuff.F90"
#include "igrid_files/timestepping_stuff.F90"
#include "igrid_files/prep_wrapup_stuff.F90"
#include "igrid_files/budgets_stuff.F90"
#include "igrid_files/popRHS_stuff.F90"
#include "igrid_files/RK45_staging.F90"

    subroutine init(this,inputfile, initialize2decomp)
        class(igrid), intent(inout), target :: this        
        character(len=clen), intent(in) :: inputfile 
        character(len=clen) :: outputdir, inputdir, scalar_info_dir, turbInfoDir, ksOutputDir, controlDir = "null", moisture_info_dir, inputDirDyaw
        integer :: nx, ny, nz, prow = 0, pcol = 0, ioUnit, nsteps = 999999, sponge_type = 1
        integer :: tid_StatsDump =10000, tid_compStats = 10000,  WallMType = 0, t_planeDump = 1000
        integer :: t_pointProbe = 10000, t_start_pointProbe = 10000, t_stop_pointProbe = 1
        integer :: runID = 0,  t_dataDump = 99999, t_restartDump = 99999,t_stop_planeDump = 1,t_dumpKSprep = 10 
        integer :: restartFile_TID = 1, ioType = 0, restartFile_RID =1, t_start_planeDump = 1
        real(rkind) :: dt=-one,tstop=one,CFL =-one,tSimStartStats=100.d0,dpfdy=zero,dPfdz=zero,CviscDT=1.d0,deltaT_dump=1.d0
        real(rkind) :: Pr = 0.7_rkind, Re = 8000._rkind, Ro = 1000._rkind,dpFdx = zero, G_alpha = 0.d0, PrandtlFluid = 1.d0, moistureFactor = 0.61_rkind
        real(rkind) :: SpongeTscale = 50._rkind, zstSponge = 0.8_rkind, Fr = 1000.d0, G_geostrophic = 1.d0
        logical ::useRestartFile=.false.,isInviscid=.false.,useCoriolis = .true., PreProcessForKS = .false.  
        logical ::isStratified=.false.,useMoisture=.false.,dumpPlanes = .false.,useExtraForcing = .false.
        logical ::useSGS = .false.,useSpongeLayer=.false.,useWindTurbines = .false., useTopAndBottomSymmetricSponge = .false. 
        logical :: useGeostrophicForcing = .false., PeriodicInZ = .false., deleteInstructions = .true., donot_dealias = .false.   
        real(rkind), dimension(:,:,:), pointer :: zinZ, zinY, zEinY, zEinZ
        integer :: AdvectionTerm = 1, NumericalSchemeVert = 0, t_DivergenceCheck = 10, ksRunID = 10
        integer :: timeSteppingScheme = 0, num_turbines = 0, P_dumpFreq = 10, P_compFreq = 10, BuoyancyTermType = 1
        logical :: normStatsByUstar=.false., ComputeStokesPressure = .true., UseDealiasFilterVert = .false., ComputeRapidSlowPressure = .false.
        real(rkind) :: tmpmn, Lz = 1.d0, latitude = 90._rkind, KSFilFact = 4.d0, dealiasFact = 2.d0/3.d0, frameAngle = 0.d0, BulkRichardson = 0.d0, HITForceTimeScale = 10.d0
        logical :: ADM = .false., storePressure = .false., useSystemInteractions = .true., &
          useFringe = .false., useHITForcing = .false., useControl = .false., &
          useHITRealSpaceLinearForcing = .false.
        integer :: tSystemInteractions = 100, ierr, KSinitType = 0, nKSvertFilt = 1, ADM_Type = 1
        logical :: computeSpectra = .false., timeAvgFullFields = .false., fastCalcPressure = .true., usedoublefringex = .false.  
        logical :: assume_fplane = .true., periodicbcs(3), useProbes = .false., KSdoZfilter = .true., computeVorticity = .false.  
        real(rkind), dimension(:,:), allocatable :: probe_locs
        real(rkind), dimension(:), allocatable :: temp
        integer :: ii, idx, temploc(1)
        logical, intent(in), optional :: initialize2decomp
        integer :: num_scalars = 0
        logical :: reset2decomp, InitSpinUp = .false., useExhaustiveFFT = .true., computeFringePressure = .false. , computeDNSPressure = .false.  
        logical :: sgsmod_stratified, dump_NU_SGS = .false., dump_KAPPA_SGS = .false., computeTurbinePressure = .false., useScalars = .false. 
        integer :: zHubIndex = 16
        real(rkind) :: angleTrigger=0.1d0, Ra = 1.d14
        character(len=4) :: scheme_xy = "FOUR"
        integer :: MeanTIDX, MeanRID, vizDump_schedule = 0    
        character(len=clen) :: MeanFilesDir, powerDumpDir 
        logical :: WriteTurbineForce = .false., useforcedStratification = .false., useDynamicYaw = .FALSE. 
        integer :: buoyancyDirection = 3, yawUpdateInterval = 100000, dealiasType = 0

        real(rkind), dimension(:,:,:), allocatable, target :: tmpzE, tmpzC, tmpyE, tmpyC
        namelist /INPUT/ nx, ny, nz, tstop, dt, CFL, nsteps, inputdir, outputdir, prow, pcol, &
                        useRestartFile, restartFile_TID, restartFile_RID, CviscDT
        namelist /IO/ vizDump_Schedule, deltaT_dump, t_restartDump, t_dataDump, ioType, dumpPlanes, runID, useProbes, &
                    & dump_NU_SGS, dump_KAPPA_SGS, t_planeDump, t_stop_planeDump, t_start_planeDump, t_start_pointProbe,&
                    & t_stop_pointProbe, t_pointProbe
        !namelist /STATS/tid_StatsDump,tid_compStats,tSimStartStats,normStatsByUstar,computeSpectra,timeAvgFullFields, computeVorticity
        namelist /PHYSICS/isInviscid,useCoriolis,useExtraForcing,isStratified,&
          useMoisture,Re,Ro,Pr,Fr, Ra, useSGS, PrandtlFluid, BulkRichardson, &
          BuoyancyTermType,useforcedStratification, useGeostrophicForcing, &
          G_geostrophic, G_alpha, dpFdx,dpFdy,dpFdz,assume_fplane,latitude,&
          useHITForcing, useScalars, frameAngle, buoyancyDirection, &
          useHITRealSpaceLinearForcing, HITForceTimeScale
        namelist /BCs/ PeriodicInZ, topWall, botWall, useSpongeLayer, zstSponge, SpongeTScale, sponge_type, botBC_Temp, topBC_Temp, useTopAndBottomSymmetricSponge, useFringe, usedoublefringex, useControl
        namelist /WINDTURBINES/ useWindTurbines, num_turbines, ADM, turbInfoDir, ADM_Type, powerDumpDir, useDynamicYaw, &
                                yawUpdateInterval, inputDirDyaw 
        namelist /NUMERICS/ AdvectionTerm, ComputeStokesPressure, NumericalSchemeVert, &
                            UseDealiasFilterVert, t_DivergenceCheck, TimeSteppingScheme, InitSpinUp, &
                            useExhaustiveFFT, dealiasFact, scheme_xy, donot_dealias, dealiasType 
        namelist /KSPREPROCESS/ PreprocessForKS, KSoutputDir, KSRunID, t_dumpKSprep, KSinitType, KSFilFact, &
                                 KSdoZfilter, nKSvertFilt
        namelist /PRESSURE_CALC/ fastCalcPressure, storePressure, P_dumpFreq, P_compFreq, computeDNSPressure, computeTurbinePressure, computeFringePressure, ComputeRapidSlowPressure            
        namelist /OS_INTERACTIONS/ useSystemInteractions, tSystemInteractions, controlDir, deleteInstructions
        namelist /SCALARS/ num_scalars, scalar_info_dir
        namelist /TURB_PRESSURE/ MeanTIDX, MeanRID, MeanFilesDir
        namelist /MOISTURE/ moistureFactor, moisture_info_dir

        ! STEP 1: READ INPUT 
        ioUnit = 11
        open(unit=ioUnit, file=trim(inputfile), form='FORMATTED', iostat=ierr)
        read(unit=ioUnit, NML=INPUT)
        read(unit=ioUnit, NML=NUMERICS)
        read(unit=ioUnit, NML=IO)
        !read(unit=ioUnit, NML=STATS)
        read(unit=ioUnit, NML=OS_INTERACTIONS)
        read(unit=ioUnit, NML=PHYSICS)
        read(unit=ioUnit, NML=PRESSURE_CALC)
        read(unit=ioUnit, NML=BCs)
        read(unit=ioUnit, NML=WINDTURBINES)
        read(unit=ioUnit, NML=KSPREPROCESS)
        this%useMoisture = useMoisture
        if (this%useMoisture) then
         read(unit=ioUnit, NML=MOISTURE)
        end if
        this%useScalars = useScalars
        if (this%useScalars) then
         read(unit=ioUnit, NML=SCALARS)
        end if
        close(ioUnit)
        this%nx = nx; this%ny = ny; this%nz = nz; this%meanfact = one/(real(nx,rkind)*real(ny,rkind)); 
        this%dt = dt; this%dtby2 = dt/two ; this%Re = Re; this%useSponge = useSpongeLayer
        this%outputdir = outputdir; this%inputdir = inputdir; this%isStratified = isStratified
        this%WallMtype = WallMType; this%runID = runID; this%tstop = tstop; this%t_dataDump = t_dataDump
        this%CFL = CFL; this%dumpPlanes = dumpPlanes; this%useGeostrophicForcing = useGeostrophicForcing
        this%timeSteppingScheme = timeSteppingScheme; this%useSystemInteractions = useSystemInteractions
        this%tSystemInteractions = tSystemInteractions; this%storePressure = storePressure
        this%P_dumpFreq = P_dumpFreq; this%P_compFreq = P_compFreq; this%timeAvgFullFields = timeAvgFullFields
        this%computeSpectra = computeSpectra; this%botBC_Temp = botBC_Temp; this%isInviscid = isInviscid
        this%assume_fplane = assume_fplane; this%useProbes = useProbes; this%PrandtlFluid = PrandtlFLuid
        this%KSinitType = KSinitType; this%KSFilFact = KSFilFact;this%useFringe = useFringe; this%useControl = useControl
        this%nsteps = nsteps; this%PeriodicinZ = periodicInZ; this%usedoublefringex = usedoublefringex 
        this%useHITForcing = useHITForcing; this%BuoyancyTermType = BuoyancyTermType; this%CviscDT = CviscDT 
        this%frameAngle = frameAngle; this%computeVorticity = computeVorticity 
        this%deleteInstructions = deleteInstructions; this%TopBC_Temp = TopBC_temp
        this%dump_NU_SGS = dump_NU_SGS; this%dump_KAPPA_SGS = dump_KAPPA_SGS; this%n_scalars = num_scalars
        this%donot_dealias = donot_dealias; this%ioType = ioType; this%HITForceTimeScale = HITForceTimeScale
        this%moistureFactor = moistureFactor; this%useHITRealSpaceLinearForcing = useHITRealSpaceLinearForcing

        if (this%CFL > zero) this%useCFL = .true. 
        if ((this%CFL < zero) .and. (this%dt < zero)) then
            call GracefulExit("Both CFL and dt cannot be negative. Have you &
            & specified either one of these in the input file?", 124)
        end if 
        this%t_restartDump = t_restartDump; this%tid_statsDump = tid_statsDump; this%useCoriolis = useCoriolis; 
        this%tSimStartStats = tSimStartStats; this%useWindTurbines = useWindTurbines
        this%tid_compStats = tid_compStats; this%useExtraForcing = useExtraForcing; this%useSGS = useSGS 
        this%UseDealiasFilterVert = UseDealiasFilterVert
        this%G_geostrophic = G_geostrophic; this%G_alpha = G_alpha; this%Fr = Fr; 
        this%fastCalcPressure = fastCalcPressure 
        this%t_start_planeDump = t_start_planeDump; this%t_stop_planeDump = t_stop_planeDump
        this%t_planeDump = t_planeDump; this%BotBC_temp = BotBC_temp; this%Ro = Ro; 
        this%PreProcessForKS = preprocessForKS; this%KSOutputDir = KSoutputDir;this%t_dumpKSprep = t_dumpKSprep 
        this%normByustar = normStatsByUstar; this%t_DivergenceCheck = t_DivergenceCheck
        this%t_start_pointProbe = t_start_pointProbe; this%t_stop_pointProbe = t_stop_pointProbe; 
        this%t_pointProbe = t_pointProbe; this%dPfdx = dPfdx; this%dPfdy = dPfdy; this%dPfdz = dPfdz
        this%InitSpinUp = InitSpinUp; this%BulkRichardson = BulkRichardson
        this%computeDNSpressure = computeDNSpressure; this%computefringePressure = computeFringePressure
        this%zHubIndex = zHubIndex; this%angleTrigger = angleTrigger
        this%computeTurbinePressure = computeTurbinePressure; this%turbPr = Pr
        this%restartPhi = 0.d0
        this%Ra = Ra
        if (useWindturbines) this%WriteTurbineForce = WriteTurbineForce

        ! STEP 2: ALLOCATE DECOMPOSITIONS
        allocate(this%gpC); allocate(this%gpE)
        if (present(initialize2decomp)) then
            reset2decomp = initialize2decomp
         else
            reset2decomp = .true.
        end if

        if (reset2decomp) then
            periodicbcs(1) = .true.; periodicbcs(2) = .true.; periodicbcs(3) = PeriodicInZ   
            call decomp_2d_init(nx, ny, nz, prow, pcol, periodicbcs)
            call get_decomp_info(this%gpC)
        else
            call decomp_info_init(nx, ny, nz, this%gpC)    
        end if

       call decomp_info_init(nx,ny,nz+1,this%gpE)
       
       if (this%useSystemInteractions) then
           if ((trim(controlDir) .eq. "null") .or.(trim(ControlDir) .eq. "NULL")) then
               this%controlDir = this%outputDir
               call message(0,"WARNING: No directory specified for OS_CONTROL instructions. Default is set to OUTPUT Directory")
           else
               this%controlDir = controlDir
           end if
       end if 
       
       if (mod(nx,2) .ne. 0) then
           call GracefulExit("The code hasn't been tested for odd values of Nx. Crazy shit could happen.", 423)
       end if 
       if (mod(ny,2) .ne. 0) then
           call GracefulExit("The code hasn't been tested for odd values of Ny. Crazy shit could happen.", 423)
       end if 
       if (mod(nz,2) .ne. 0) then
           call GracefulExit("The code hasn't been tested for odd values of Nz. Crazy shit could happen.", 423)
       end if 

       call get_boundary_conditions_stencil()

       ! Set numerics
       select case(NumericalSchemeVert)
       case(0)
           useCompactFD = .false.
       case(1)
           useCompactFD = .true.
       case(2)
           useCompactFD = .false.
           if (.not. PeriodicInZ) then
              call gracefulExit("If you use Fourier Collocation in Z, the problem must be periodic in Z.",123)
           end if
       case default
           call gracefulExit("Invalid choice for NUMERICALSCHEMEVERT",423)
       end select

       select case(AdvectionTerm)
       case(0)
           useSkewSymm = .false.
       case(1)
           useSkewSymm = .true.
       case default
           call gracefulExit("Invalid choice for ADVECTIONTERM",423)
       end select
       
       ! STEP 3: GENERATE MESH (CELL CENTERED) 
       if ( allocated(this%mesh) ) deallocate(this%mesh) 
       allocate(this%mesh(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3),3))
       call meshgen_WallM(this%gpC, this%dx, this%dy, &
           this%dz, this%mesh,inputfile) ! <-- this procedure is part of user defined HOOKS
       this%zTop = p_maxval(this%mesh(:,:,:,3)) + this%dz/2.d0
       this%zBot = p_minval(this%mesh(:,:,:,3)) - this%dz/2.d0
       this%zMid = half*(this%zTop + this%zBot)
       Lz = this%zTop - this%zBot
       call message(0,"Mesh generated:")
       call message(1,"dx:", this%dx)
       call message(1,"dy:", this%dy)
       call message(1,"dz:", this%dz)
       call message(1,"Lz:", Lz)


       ! STEP 4: ALLOCATE/INITIALIZE THE SPECTRAL DERIVED TYPES
       allocate(this%spectC)
       call this%spectC%init("x", nx, ny, nz  , this%dx, this%dy, this%dz, &
               scheme_xy, this%filter_x, 2 , fixOddball=.false., exhaustiveFFT=useExhaustiveFFT, init_periodicInZ=periodicinZ, dealiasF=dealiasfact, dealiasType=dealiasType)
       allocate(this%spectE)
       call this%spectE%init("x", nx, ny, nz+1, this%dx, this%dy, this%dz, &
               scheme_xy, this%filter_x, 2 , fixOddball=.false., exhaustiveFFT=useExhaustiveFFT, init_periodicInZ=.false., dealiasF=dealiasfact, dealiasType=dealiasType)
       this%sp_gpC => this%spectC%spectdecomp
       this%sp_gpE => this%spectE%spectdecomp


       ! STEP 5: ALLOCATE/INITIALIZE THE DERIVATIVE DERIVED TYPE
       allocate(this%Pade6OpZ)
       call this%Pade6OpZ%init(this%gpC,this%sp_gpC, this%gpE, this%sp_gpE,this%dz,NumericalSchemeVert,PeriodicInZ,this%spectC)
       allocate(this%OpsPP)
       call this%OpsPP%init(this%gpC,this%gpE,0,this%dx,this%dy,this%dz,this%spectC%spectdecomp, &
                   this%spectE%spectdecomp, .false., .false.)
       
       if (this%UseDealiasFilterVert) then
           allocate(this%filzC, this%filzE)
           ierr = this%filzC%init(nz  , PeriodicInZ)
           ierr = this%filzE%init(nz+1, PeriodicInZ)
       end if


       ! STEP 6: ALLOCATE MEMORY FOR FIELD ARRAYS
       allocate(this%PfieldsC(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3),7))
       allocate(this%PfieldsE(this%gpE%xsz(1),this%gpE%xsz(2),this%gpE%xsz(3),4))
       allocate(this%dTdzE(this%gpE%xsz(1),this%gpE%xsz(2),this%gpE%xsz(3)))
       allocate(this%dTdzC(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)))
       allocate(this%dTdxC(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)))
       allocate(this%dTdyC(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)))
       allocate(this%dTdxE(this%gpE%xsz(1),this%gpE%xsz(2),this%gpE%xsz(3)))
       allocate(this%dTdyE(this%gpE%xsz(1),this%gpE%xsz(2),this%gpE%xsz(3)))
       call this%spectC%alloc_r2c_out(this%SfieldsC,4)
       call this%spectC%alloc_r2c_out(this%dTdxH)
       call this%spectC%alloc_r2c_out(this%dTdyH)
       call this%spectE%alloc_r2c_out(this%dTdzH)
       call this%spectC%alloc_r2c_out(this%dTdzHC)
       call this%spectC%alloc_r2c_out(this%rhsC,3) 
       call this%spectC%alloc_r2c_out(this%OrhsC,3)
       call this%spectE%alloc_r2c_out(this%SfieldsE,4)
       allocate(this%divergence(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)))
       allocate(this%duidxjC(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3),9))
       allocate(this%duidxjE(this%gpE%xsz(1),this%gpE%xsz(2),this%gpE%xsz(3),9))
       call this%spectC%alloc_r2c_out(this%duidxjChat,9)
       call this%spectE%alloc_r2c_out(this%duidxjEhat,9)
       call this%spectE%alloc_r2c_out(this%rhsE,1); call this%spectE%alloc_r2c_out(this%OrhsE,1)
       
       this%u => this%PfieldsC(:,:,:,1) ; this%v => this%PfieldsC(:,:,:,2) ; this%wC => this%PfieldsC(:,:,:,3) 
       this%w => this%PfieldsE(:,:,:,1) ; this%uE => this%PfieldsE(:,:,:,2) ; this%vE => this%PfieldsE(:,:,:,3) 
       
       this%uhat => this%SfieldsC(:,:,:,1); this%vhat => this%SfieldsC(:,:,:,2); 
       this%whatC => this%SfieldsC(:,:,:,3); this%what => this%SfieldsE(:,:,:,1)
       this%uEhat => this%SfieldsE(:,:,:,3); this%vEhat => this%SfieldsE(:,:,:,4)

       this%ox => this%PfieldsC(:,:,:,4); this%oy => this%PfieldsC(:,:,:,5); this%oz => this%PfieldsC(:,:,:,6)

       this%u_rhs => this%rhsC(:,:,:,1); this%v_rhs => this%rhsC(:,:,:,2); this%w_rhs => this%rhsE(:,:,:,1)

       this%u_Orhs => this%OrhsC(:,:,:,1); this%v_Orhs => this%OrhsC(:,:,:,2); this%w_Orhs => this%OrhsE(:,:,:,1)

       !if (this%isStratified) then
       this%T => this%PfieldsC(:,:,:,7); this%That => this%SfieldsC(:,:,:,4)
       this%TE => this%PfieldsE(:,:,:,4); this%T_rhs => this%rhsC(:,:,:,3)
       this%T_Orhs => this%OrhsC(:,:,:,3); this%TEhat => this%SfieldsE(:,:,:,2)
       !end if

       !allocate(this%cbuffxC(this%sp_gpC%xsz(1),this%sp_gpC%xsz(2),this%sp_gpC%xsz(3),2))
       allocate(this%cbuffyC(this%sp_gpC%ysz(1),this%sp_gpC%ysz(2),this%sp_gpC%ysz(3),3))
       allocate(this%cbuffyE(this%sp_gpE%ysz(1),this%sp_gpE%ysz(2),this%sp_gpE%ysz(3),2))
       
       allocate(this%cbuffzC(this%sp_gpC%zsz(1),this%sp_gpC%zsz(2),this%sp_gpC%zsz(3),3))
       allocate(this%cbuffzE(this%sp_gpE%zsz(1),this%sp_gpE%zsz(2),this%sp_gpE%zsz(3),2))

       allocate(this%rbuffxC(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3),4))
       allocate(this%rbuffyC(this%gpC%ysz(1),this%gpC%ysz(2),this%gpC%ysz(3),2))
       allocate(this%rbuffzC(this%gpC%zsz(1),this%gpC%zsz(2),this%gpC%zsz(3),4))

       allocate(this%rbuffxE(this%gpE%xsz(1),this%gpE%xsz(2),this%gpE%xsz(3),2))
       allocate(this%rbuffyE(this%gpE%ysz(1),this%gpE%ysz(2),this%gpE%ysz(3),2))
       allocate(this%rbuffzE(this%gpE%zsz(1),this%gpE%zsz(2),this%gpE%zsz(3),4))
       if (.not.this%isinviscid) then
           allocate(this%d2udz2hatC(this%sp_gpC%ysz(1),this%sp_gpC%ysz(2),this%sp_gpC%ysz(3)))
           allocate(this%d2vdz2hatC(this%sp_gpC%ysz(1),this%sp_gpC%ysz(2),this%sp_gpC%ysz(3)))
           allocate(this%d2wdz2hatE(this%sp_gpE%ysz(1),this%sp_gpE%ysz(2),this%sp_gpE%ysz(3)))
           if (this%isStratified) then
              allocate(this%d2Tdz2hatC(this%sp_gpC%ysz(1),this%sp_gpC%ysz(2),this%sp_gpC%ysz(3)))
           end if 
       end if 
       this%nxZ = size(this%cbuffzE,1); this%nyZ = size(this%cbuffzE,2)
       allocate(this%fbody_x(this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3)))
       allocate(this%fbody_y(this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3)))
       allocate(this%fbody_z(this%gpE%xsz(1), this%gpE%xsz(2), this%gpE%xsz(3)))
       this%storeFbody = .true. ! Cant think of a case where this will be false 


       ! STEP 6: ALLOCATE/INITIALIZE THE POISSON DERIVED TYPE
       allocate(this%padepoiss)
       call this%padepoiss%init(this%dx , this%dy, this%dz, this%spectC, this%spectE, computeStokesPressure, Lz, .true., &
                                this%gpC, this%Pade6opz, PeriodicInZ) 
          
       ! STEP 7: INITIALIZE THE FIELDS
       if (useRestartFile) then
           call this%readRestartFile(restartfile_TID, restartfile_RID)
           this%step = restartfile_TID
       else 
           call initfields_wallM(this%gpC, this%gpE, inputfile, this%mesh, this%PfieldsC, this%PfieldsE)! <-- this procedure is part of user defined HOOKS
           this%step = 0
           this%tsim = zero
           !call this%dumpRestartfile()
       end if 
   
       if (this%isStratified) then
           if (BuoyancyTermType == 1) then
               call set_Reference_Temperature(inputfile,this%ThetaRef)
               call message(1,"Reference Temperature set to:",this%ThetaRef) 
           end if
          
           if (botBC_Temp == 0) then
               call setDirichletBC_Temp(inputfile, this%Tsurf, this%dTsurf_dt)
               this%Tsurf0 = this%Tsurf
               this%Tsurf = this%Tsurf0 + this%dTsurf_dt*this%tsim
           else if (botBC_Temp == 1) then
               ! Do nothing 
           else if (botBC_Temp == 2) then
               this%wTh_surf = this%tsim
               call setInhomogeneousNeumannBC_Temp(inputfile, this%wTh_surf)
           else if (botBC_Temp == 3) then
               if (topBC_Temp .ne. 3) then
                    call GraceFulExit("Zero dirichlet BC must be set symmetrically if botBC_Temp = 3",341)        
               end if 
           else
               call GraceFulExit("Only Dirichlet, Homog. Neumann or Inhomog. Neumann BCs supported for Temperature at &
                   & this time. Set botBC_Temp = 0 or 1 or 2",341)        
           end if
           if (topBC_Temp == 0) then
               call setDirichletBC_Temp(inputfile, this%TsurfTop, this%dTsurfTop_dt)
               this%Tsurf0 = this%TsurfTop
               this%TsurfTop = this%Tsurf0 + this%dTsurfTop_dt*this%tsim
           end if 
        end if 

       if (this%initspinup) then
          if (this%isStratified) then
           call GracefulExit("InitSpinUp not permitted when stratification is ON",3124)
          else
           call this%Initialize_Scalar_for_InitSpinUp(useRestartFile, inputfile, restartfile_TID, restartfile_RID)
          end if
       end if

       call this%spectC%fft(this%u,this%uhat)   
       call this%spectC%fft(this%v,this%vhat)   
       call this%spectE%fft(this%w,this%what)   
       if (this%isStratified .or. this%initspinup) call this%spectC%fft(this%T,this%That)   

       ! Dealias and filter before projection
       call this%dealiasFields()


       ! Pressure projection
       call this%padepoiss%DivergenceCheck(this%uhat, this%vhat, this%what, this%divergence)
       call this%padepoiss%PressureProjection(this%uhat,this%vhat,this%what)
       !call this%padepoiss%DivergenceCheck(this%uhat, this%vhat, this%what, this%divergence,.true.)

       ! Take it back to physical fields
       call this%spectC%ifft(this%uhat,this%u)
       call this%spectC%ifft(this%vhat,this%v)
       call this%spectE%ifft(this%what,this%w)
       if (this%isStratified) call this%spectC%ifft(this%That,this%T)

       ! STEP 8: Interpolate the cell center values of w
       !if (this%useSGS) then
       !    call this%compute_and_bcast_surface_Mn()
       !end if

       !if ((PeriodicInZ) .and. (useHITforcing)) then
       !    tmpmn = p_sum(this%u)/(real(nx,rkind)*real(ny,rkind)*real(nz,rkind))
       !    this%u = this%u - tmpmn
       !    
       !    tmpmn = p_sum(this%v)/(real(nx,rkind)*real(ny,rkind)*real(nz,rkind))
       !    this%v = this%v - tmpmn
       !    
       !    tmpmn = p_sum(this%w)/(real(nx,rkind)*real(ny,rkind)*real(nz + 1,rkind))
       !    this%w = this%w - tmpmn
       !end if
       call this%interp_PrimitiveVars()
       call message(1,"Max KE:",P_MAXVAL(this%getMaxKE()))
    
       ! STEP 9: Compute duidxj
       call this%compute_duidxj()
       if (this%isStratified) call this%compute_dTdxi() 

       ! STEP 10a: Compute Coriolis Term
       if (this%useCoriolis) then
           call message(0, "Turning on Coriolis with Geostrophic Forcing")
           call message(1, "Geostrophic Velocity Magnitude  :", this%G_geostrophic) 
           call message(1, "Geostrophic Velocity Direction  :", this%G_alpha)
           call message(1, "Rossby Number:", this%Ro) 
           call this%spectC%alloc_r2c_out(this%GxHat)
           call this%spectC%alloc_r2c_out(this%GyHat)
           call this%spectE%alloc_r2c_out(this%GxHat_Edge)
           call this%spectE%alloc_r2c_out(this%GyHat_Edge)
           this%rbuffxC(:,:,:,1) = this%G_GEOSTROPHIC*cos(G_ALPHA*pi/180.d0)
           call this%spectC%fft(this%rbuffxC(:,:,:,1),this%Gxhat)
           this%rbuffxE(:,:,:,1) = this%G_GEOSTROPHIC*cos(G_ALPHA*pi/180.d0)
           call this%spectE%fft(this%rbuffxE(:,:,:,1),this%Gxhat_Edge)
           this%rbuffxC(:,:,:,1) = this%G_GEOSTROPHIC*sin(G_ALPHA*pi/180.d0)
           call this%spectC%fft(this%rbuffxC(:,:,:,1),this%Gyhat)
           this%rbuffxE(:,:,:,1) = this%G_GEOSTROPHIC*sin(G_ALPHA*pi/180.d0)
           call this%spectE%fft(this%rbuffxE(:,:,:,1),this%Gyhat_Edge)

           this%latitude = latitude
           this%frameAngle = frameAngle
           
           if (this%assume_fplane) then
               this%coriolis_omegaZ   = sin(latitude*pi/180.d0)
               this%coriolis_omegaY = 0.d0
               this%coriolis_omegaX = 0.d0
               call message(1, "Making the f-plane assumption (latitude effect &
               & ignored in w equation)")
           else
               this%coriolis_omegaX = cos(latitude*pi/180.d0)*sin(frameAngle*pi/180.d0)
               this%coriolis_omegaZ = sin(latitude*pi/180.d0)
               this%coriolis_omegaY = cos(latitude*pi/180.d0)*cos(frameAngle*pi/180.d0)
               call message(1,"Latitude used for Coriolis (degrees)",latitude)
           end if
       end if

       ! STEP 10b: Compute additional forcing (channel)
       if (this%useExtraForcing) then
           call message(0," Turning on aditional forcing")
           call message(1," dP_dx = ", dpFdx)
           call this%spectC%alloc_r2c_out(this%dpF_dxhat)
           this%rbuffxC(:,:,:,1) = dpFdx
           call this%spectC%fft(this%rbuffxC(:,:,:,1),this%dpF_dxhat)
       end if  

        ! STEP 11: Initialize SGS model
        allocate(this%SGSmodel)
        if (this%useSGS) then
            ! First get z at edges
            zinY => this%rbuffyC(:,:,:,1); zinZ => this%rbuffzC(:,:,:,1)
            zEinZ => this%rbuffzE(:,:,:,1); zEinY => this%rbuffyE(:,:,:,1)
            call transpose_x_to_y(this%mesh(:,:,:,3),zinY,this%gpC)
            call transpose_y_to_z(zinY,zinZ,this%gpC)
            call this%OpsPP%InterpZ_Cell2Edge(zinZ,zEinZ,zero,zero)
            zEinZ(:,:,this%nz+1) = zEinZ(:,:,this%nz) + this%dz
            call transpose_z_to_y(zEinZ,zEinY,this%gpE)
            call transpose_y_to_x(zEinY,this%rbuffxE(:,:,:,1), this%gpE)

            if ((this%initSpinup) .or. (this%useScalars) .or. (this%isStratified)) then
               sgsmod_stratified = .true. 
            else
               sgsmod_stratified = .false. 
            end if 
            call this%sgsModel%init(this%gpC, this%gpE, this%spectC, this%spectE, this%dx, this%dy, this%dz, inputfile, &
                                    this%rbuffxE(1,1,:,1), this%mesh(1,1,:,3), this%fBody_x, this%fBody_y, this%fBody_z, &
                                    this%storeFbody,this%Pade6opZ, this%cbuffyC, this%cbuffzC, this%cbuffyE, this%cbuffzE, &
                                    this%rbuffxC, this%rbuffyC, this%rbuffzC, this%rbuffyE, this%rbuffzE, this%Tsurf, &
                                    this%ThetaRef, this%wTh_surf, this%Fr, this%Re, this%isInviscid, sgsmod_stratified, &
                                    this%botBC_Temp, this%initSpinUp)
            call this%sgsModel%link_pointers(this%nu_SGS, this%tauSGS_ij, this%tau13, this%tau23, this%q1, this%q2, this%q3, this%kappaSGS)
            call message(0,"SGS model initialized successfully")
        end if 
        this%max_nuSGS = zero


        ! STEP 12: Set Sponge Layer
        if (this%useSponge) then
            allocate(this%RdampC(this%sp_gpC%ysz(1), this%sp_gpC%ysz(2), this%sp_gpC%ysz(3)))
            allocate(this%RdampE(this%sp_gpE%ysz(1), this%sp_gpE%ysz(2), this%sp_gpE%ysz(3)))
            zinZ => this%rbuffzC(:,:,:,1)
            zEinZ => this%rbuffzE(:,:,:,1); 
            call transpose_x_to_y(this%mesh(:,:,:,3),zinY,this%gpC)
            call transpose_y_to_z(zinY,zinZ,this%gpC)
            call this%OpsPP%InterpZ_Cell2Edge(zinZ,zEinZ,zero,zero)
            zEinZ(:,:,this%nz+1) = zEinZ(:,:,this%nz) + this%dz
            !ztop = zEinZ(1,1,this%nz+1); 
            if (zstSponge >= 1) then
                call GracefulExit("zstSponge must be less than 1.",245)
            end if
            
            allocate(tmpyC(this%sp_gpC%ysz(1),this%sp_gpC%ysz(2),this%sp_gpC%ysz(3)))
            allocate(tmpyE(this%sp_gpE%ysz(1),this%sp_gpE%ysz(2),this%sp_gpE%ysz(3)))

            zinY => tmpyC
            zEinY => tmpyE
            allocate(tmpzC(this%sp_gpC%zsz(1),this%sp_gpC%zsz(2),this%sp_gpC%zsz(3)))
            do idx = 1,size(zinZ,3)
                tmpzC(:,:,idx) = zinZ(1,1,idx)
            end do 
            call transpose_z_to_y(tmpzC,zinY,this%sp_gpC)
            deallocate(tmpzC)
            allocate(tmpzE(this%sp_gpE%zsz(1),this%sp_gpE%zsz(2),this%sp_gpE%zsz(3)))
            do idx = 1,size(zEinZ,3)
                tmpzE(:,:,idx) = zEinZ(1,1,idx)
            end do 
            call transpose_z_to_y(tmpzE,zEinY, this%sp_gpC) 
            deallocate(tmpzE)
            nullify(zEinZ, zinZ)

            zstSponge = zstSponge*(this%zTop - this%zBot)     !! <PERCENTAGE OF THE DOMAIN>
            !zstSponge = zstSponge*this%zTop                   !! <PERCENTAGE OF THE DOMAIN>
            select case(sponge_type)
            case(1)

                !!!!!!!!!!!!!!! OLD BLOCK !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !this%RdampC = (one/SpongeTscale) * (one - cos(pi*(zinY - zstSponge) /(this%zTop - zstSponge)))/two
                !this%RdampE = (one/SpongeTscale) * (one - cos(pi*(zEinY - zstSponge)/(this%zTop - zstSponge)))/two
                !if (useTopAndBottomSymmetricSponge) then
                !   where (abs(zEinY) < zstSponge) 
                !       this%RdampE = zero
                !   end where
                !   where (abs(zinY) < zstSponge) 
                !       this%RdampC = zero
                !   end where
                !   if ((zEinZ(1,1,1) + zEinZ(1,1,this%nz+1))<1.d-13) then
                !      call message(0,"WARNING: Computed domain is not symmetric &
                !         & about z=0. You shouldn't use the symmetric sponge")
                !      call MPI_BARRIER(mpi_comm_world, ierr)
                !      call GracefulExit("Failed at sponge initialization",134)
                !   end if   
                !else
                !   where (zEinY < zstSponge) 
                !       this%RdampE = zero
                !   end where
                !   where (zinY < zstSponge) 
                !       this%RdampC = zero
                !   end where
                !end if 
                !!!!!!!!!!!!!!! OLD BLOCK -- assumes zBot = 0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                
                if (useTopAndBottomSymmetricSponge) then
                    ! Ensure zinY and zEinY are centered at 0 irrespective of this%zBot
                    zinY  = zinY  - this%zMid
                    zEinY = zEinY - this%zMid

                    this%RdampC = (one/SpongeTscale) * (one - cos(pi*(abs(zinY) - zstSponge + this%zMid) /(this%zTop - zstSponge)))/two
                    this%RdampE = (one/SpongeTscale) * (one - cos(pi*(abs(zEinY) - zstSponge + this%zMid)/(this%zTop - zstSponge)))/two

                    where (abs(zEinY) < (zstSponge-this%zMid)) 
                        this%RdampE = zero
                    end where
                    where (abs(zinY) < (zstSponge-this%zMid)) 
                        this%RdampC = zero
                    end where
                    !call this%dump(this%RdampC,'spng',this%gpC)
                    !call this%dump(this%RdampE,'spng',this%gpE)
                else
                    ! Ensure zinY and zEinY start at 0 irrespective of this%zBot
                    zinY  = zinY  - this%zBot
                    zEinY = zEinY - this%zBot

                    this%RdampC = (one/SpongeTscale) * (one - cos(pi*(zinY - zstSponge) /(this%zTop - zstSponge)))/two
                    this%RdampE = (one/SpongeTscale) * (one - cos(pi*(zEinY - zstSponge)/(this%zTop - zstSponge)))/two

                   where (zEinY < zstSponge) 
                       this%RdampE = zero
                   end where
                   where (zinY < zstSponge) 
                       this%RdampC = zero
                   end where
                end if 
            case(2)
                zinY = (zinY - zstSponge)/(2*(this%zTop - zstSponge))
                call S_sponge_smooth(zinY,this%RdampC)
                this%RdampC = (2.5d0/SpongeTscale)*this%RdampC
                
                zEinY = (zEinY - zstSponge)/(2*(this%zTop - zstSponge))
                call S_sponge_smooth(zEinY,this%RdampE)
                this%RdampE = (2.5d0/SpongeTscale)*this%RdampE ! consistent noramlization with type 1 sponge
            end select 

            nullify(zinY, zEinY)
            deallocate(tmpyC, tmpyE)
            ! REMOVING THE LINES BELOW TO ENSURE THAT SPONGE ONLY AFFECTS
            ! FLUCTUATIONS
            !call this%spectC%alloc_r2c_out(this%uBase)
            !call this%spectC%alloc_r2c_out(this%vBase)
            !call this%spectC%alloc_r2c_out(this%TBase)
            !this%rbuffxC(:,:,:,1) = this%u 
            !call this%spectC%fft(this%rbuffxC(:,:,:,1),this%uBase)
            !this%rbuffxC(:,:,:,1) = this%v
            !call this%spectC%fft(this%rbuffxC(:,:,:,1),this%vBase)
            !this%rbuffxC(:,:,:,1) = this%T
            !call this%spectC%fft(this%rbuffxC(:,:,:,1),this%TBase)
            call message(0,"Sponge Layer initialized successfully")
            call message(1,"Sponge Layer active above z = ",zstSponge)
        end if 
       
        if (this%useWindTurbines) then
           allocate(this%WindTurbineArr)
           call this%WindTurbineArr%init(inputFile, this%gpC, this%gpE, this%spectC, this%spectE, this%cbuffyC, this%cbuffyE, this%cbuffzC, this%cbuffzE, this%mesh, this%dx, this%dy, this%dz)
           allocate(this%inst_horz_avg_turb(8*this%WindTurbineArr%nTurbines))
           this%inst_horz_avg_turb = zero
       end if
       allocate(this%inst_horz_avg(5))
       this%inst_horz_avg = zero

       ! STEP 13: Set visualization planes for io
       call set_planes_io(this%xplanes, this%yplanes, this%zplanes)


       ! STEP 14a : Probes
       if (this%useProbes) then
           call hook_probes(inputfile, probe_locs)
           if (.not. allocated(probe_locs)) then
               call GracefulExit("You forgot to set the probe locations in initialize.F90 file for the problem.",123)
           end if
           if (size(probe_locs,2) > 999) then
               call GracefulExit("Maximum number of probes allowed is 999",123)
           end if

           this%nprobes = 0
           do idx = 1,size(probe_locs,2)
               ! assume x - decomposition
               ! First check if y lies within decomposition
               if ((probe_locs(2,idx) < maxval(this%mesh(:,:,:,2))+this%dy/2.d0) .and. (probe_locs(2,idx) > minval(this%mesh(:,:,:,2))- this%dy/2.d0)) then
                   ! Now check if z lies within my decomposition
                   if ((probe_locs(3,idx) < maxval(this%mesh(:,:,:,3))+this%dz/2.d0) .and. (probe_locs(3,idx) > minval(this%mesh(:,:,:,3))-this%dz/2.d0)) then
                       ! Looks like I have the probe!
                       this%nprobes = this%nprobes + 1
                   end if
               end if
           end do  
              
           ! If have 1 or more probes, I need to allocate memory for probes
           if (this%nprobes > 0) then
               allocate(this%probes(4,this%nprobes))
               !if (this%isStratified) then
               !    allocate(this%probe_data(1:5,1:this%nprobes,0:this%probeTimeLimit-1)) ! Store time + 3 fields
               !else                                     
               !    allocate(this%probe_data(1:4,1:this%nprobes,0:this%probeTimeLimit-1)) ! Store time + 3 fields
               !end if
               allocate(this%probe_data(1:9,1:this%nprobes,0:this%probeTimeLimit-1))
               this%probe_data = 0.d0
               ii = 1
               do idx = 1,size(probe_locs,2)
                   if ((probe_locs(2,idx) < maxval(this%mesh(:,:,:,2)) + this%dy/2.d0) .and. (probe_locs(2,idx) > minval(this%mesh(:,:,:,2)) - this%dy/2.d0)) then
                       if ((probe_locs(3,idx) < maxval(this%mesh(:,:,:,3)) + this%dz/2.d0) .and. (probe_locs(3,idx) > minval(this%mesh(:,:,:,3)) - this%dz/2.d0)) then
                           allocate(temp(size(this%mesh,1)))
                           temp = abs(probe_locs(1,idx) - this%mesh(:,1,1,1))
                           temploc = minloc(temp)
                           this%probes(1,ii) = temploc(1) 
                           deallocate(temp)

                           allocate(temp(size(this%mesh,2)))
                           temp = abs(probe_locs(2,idx) - this%mesh(1,:,1,2));
                           temploc = minloc(temp)
                           this%probes(2,ii) = temploc(1)
                           deallocate(temp)
                           
                           allocate(temp(size(this%mesh,3)))
                           temp = abs(probe_locs(3,idx) - this%mesh(1,1,:,3))
                           temploc = minloc(temp)
                           this%probes(3,ii) = temploc(1) 
                           deallocate(temp)
                           
                           this%probes(4,ii) = idx ! Probe ID
                           ii = ii + 1
                       end if
                   end if
               end do 
               this%doIhaveAnyProbes = .true. 
           else
               this%doIhaveAnyProbes = .false.  
           end if
           this%ProbeStartStep = this%step
           deallocate(probe_locs)
           !print*, nrank, "Do I have probes?:", this%doIhaveAnyProbes, this%nprobes
           call message(0,"Total probes initialized:", p_sum(this%nprobes))
       end if
     
       ! STEP 14b : Preprocessing for KS
       if (this%PreprocessForKS) then
           allocate(this%LES2KS)
           if (this%KSinitType == 0) then
               !call this%LES2KS%init(this%spectC, this%gpC, this%dx, this%dy, this%outputdir, this%RunID, this%probes, this%KSFilFact, KSdoZfilter, nKSvertFilt)
               call this%LES2KS%init(this%spectC, this%gpC, this%dx, this%dy, this%KSOutputDir, this%RunID, this%probes, this%KSFilFact, KSdoZfilter, nKSvertFilt)
               call this%LES2KS%link_pointers(this%uFil4KS, this%vFil4KS, this%wFil4KS)
               call this%LES2KS%applyFilterForKS(this%u, this%v, this%w)
               call this%LES2KS%dumpKSfilteredFields()
               if (this%useProbes) then
                   if (this%doIhaveAnyProbes) then
                       allocate(this%KS_Probe_Data(1:4,1:this%nprobes,0:this%probeTimeLimit-1))
                   end if
               end if
           else
               call GracefulExit("All KSinitTypes except for 0 are temporarily suspended.",12)
               call set_KS_planes_io(this%planes2dumpC_KS, this%planes2dumpF_KS) 
               call this%LES2KS%init(nx,ny,nz,this%spectE, this%gpE, this%KSOutputDir, KSrunID, this%dx, this%dy, &
                  &         this%dz, this%planes2dumpC_KS, this%planes2dumpF_KS)
           end if
           this%KSupdated = .false. 
           call message(0, "KS Preprocessor initialized successfully.")
       end if 


       ! STEP 15: Set up extra buffers for RK3
       if (timeSteppingScheme == 1) then
           if (this%isStratified .or. this%initspinup) then
               call this%spectC%alloc_r2c_out(this%SfieldsC2,3)
               call this%spectE%alloc_r2c_out(this%SfieldsE2,2)
           else
               call this%spectC%alloc_r2c_out(this%SfieldsC2,2)
               call this%spectE%alloc_r2c_out(this%SfieldsE2,1)
           end if 
           this%uhat1 => this%SfieldsC2(:,:,:,1); 
           this%vhat1 => this%SfieldsC2(:,:,:,2); 
           this%what1 => this%SfieldsE2(:,:,:,1); 
           if (this%isStratified .or. this%initspinup) then
               this%That1 => this%SfieldsC2(:,:,:,3); 
           end if 
        else if (timeSteppingScheme == 2) then
           call this%spectC%alloc_r2c_out(this%uExtra,3)
           call this%spectC%alloc_r2c_out(this%uRHSExtra,1)
           call this%spectC%alloc_r2c_out(this%vExtra,3)
           call this%spectC%alloc_r2c_out(this%vRHSExtra,1)
           call this%spectE%alloc_r2c_out(this%wExtra,3)
           call this%spectE%alloc_r2c_out(this%wRHSExtra,1)
           call this%spectC%alloc_r2c_out(this%TExtra,3)
           call this%spectC%alloc_r2c_out(this%TRHSExtra,1)
        end if 

       if ((timeSteppingScheme .ne. 0) .and. (timeSteppingScheme .ne. 1) .and. (timeSteppingScheme .ne. 2)) then
           call GracefulExit("Invalid choice of TIMESTEPPINGSCHEME.",5235)
       end if 

       ! STEP 16: Initialize Statistics
       if (this%timeAvgFullFields) then
           call this%init_stats3D()
       else
       !    call this%init_stats()
       end if
      
       ! STEP 17: Set Fringe
       allocate(this%fringe_x1, this%fringe_x2)
       allocate(this%fringe_x)
       if (this%usedoublefringex) then
           call this%fringe_x1%init(inputfile, this%dx, this%mesh(:,1,1,1), this%dy, this%mesh(1,:,1,2), &
                                       this%dz, this%mesh(1,1,:,3), & 
                                       this%spectC, this%spectE, this%gpC, this%gpE, &
                                       this%rbuffxC, this%rbuffxE, this%cbuffyC, this%cbuffyE, fringeID=1)   
           call this%fringe_x2%init(inputfile, this%dx, this%mesh(:,1,1,1), this%dy, this%mesh(1,:,1,2), &
                                       this%dz, this%mesh(1,1,:,3), & 
                                       this%spectC, this%spectE, this%gpC, this%gpE, &
                                       this%rbuffxC, this%rbuffxE, this%cbuffyC, this%cbuffyE, fringeID=2)   
       else
           if (this%useFringe) then
               call this%fringe_x%init(inputfile, this%dx, this%mesh(:,1,1,1), this%dy, this%mesh(1,:,1,2), &
                                       this%dz, this%mesh(1,1,:,3), & 
                                       this%spectC, this%spectE, this%gpC, this%gpE, &
                                       this%rbuffxC, this%rbuffxE, this%cbuffyC, this%cbuffyE)   
           end if
       end if 
       
       ! STEP 18: Set HIT Forcing
       if (this%useHITForcing) then
           allocate(this%hitforce)
           call this%hitforce%init(inputfile, this%gpC, this%sp_gpC, this%sp_gpE, &
             this%spectC, this%spectE, this%cbuffyE(:,:,:,1), this%cbuffyC(:,:,:,1), &
             this%cbuffzE(:,:,:,1), this%cbuffzC, this%rbuffxC(:,:,:,1), this%step, &
             this%padepoiss)
       end if
       
       ! STEP 19: Set up storage for Pressure
       if ((this%storePressure) .or. (this%fastCalcPressure)) then
           allocate(this%Pressure(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)))
           if (this%computefringePressure) then
              allocate(this%urhs_fringe(this%sp_gpC%ysz(1),this%sp_gpC%ysz(2),this%sp_gpC%ysz(3)))
              allocate(this%vrhs_fringe(this%sp_gpC%ysz(1),this%sp_gpC%ysz(2),this%sp_gpC%ysz(3)))
              allocate(this%wrhs_fringe(this%sp_gpE%ysz(1),this%sp_gpE%ysz(2),this%sp_gpE%ysz(3)))
              allocate(this%pressure_fringe(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)))
           end if
           if (this%computeDNSpressure) then
              allocate(this%urhs_dns(this%sp_gpC%ysz(1),this%sp_gpC%ysz(2),this%sp_gpC%ysz(3)))
              allocate(this%vrhs_dns(this%sp_gpC%ysz(1),this%sp_gpC%ysz(2),this%sp_gpC%ysz(3)))
              allocate(this%wrhs_dns(this%sp_gpE%ysz(1),this%sp_gpE%ysz(2),this%sp_gpE%ysz(3)))
              allocate(this%pressure_dns(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)))
           end if
           if (this%computeTurbinePressure) then
              allocate(this%urhs_turbine(this%sp_gpC%ysz(1),this%sp_gpC%ysz(2),this%sp_gpC%ysz(3)))
              allocate(this%vrhs_turbine(this%sp_gpC%ysz(1),this%sp_gpC%ysz(2),this%sp_gpC%ysz(3)))
              allocate(this%wrhs_turbine(this%sp_gpE%ysz(1),this%sp_gpE%ysz(2),this%sp_gpE%ysz(3)))
              allocate(this%pressure_turbine(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)))
              this%urhs_turbine = dcmplx(0.d0, 0.d0)
              this%vrhs_turbine = dcmplx(0.d0, 0.d0)
              this%wrhs_turbine = dcmplx(0.d0, 0.d0)
           end if
           call message(1, "Done allocating storage for pressure")
       end if
       if ((.not. this%fastCalcPressure) .and. ((this%computefringePressure) .or. (this%computeDNSPressure))) then
           call gracefulExit("You need to set FASTCALCPRESSURE = .true. in & 
                    & order to use computefringepressure or computeDNSpressure", 313)
       end if

       ! STEP 20: Update the probes
       if (this%useProbes) call this%updateProbes()

       ! STEP 21: Buoyancy term type
       if (this%isStratified) then
           select case (this%BuoyancyTermType)
           case(1)
              ! Parameter : Froude number 
              ! RHS_3 = RHS_3 + one/(Fr^2*Theta_ref)*(T - T_ref) ! implemented with T_ref absorbed as hydrostatic pressure
              this%BuoyancyFact = one/(this%Fr*this%Fr*this%ThetaRef)
              call message(1,"Buoyancy term type 1 selected. Buoyancy term &
                                & calculation term uses")
              call message(2,"Froude number:", this%Fr)
              call message(2,"Reference temperature:", this%thetaRef)
           case(2)
                ! Parameter : Bulk Richardson number 
                ! RHS_3 = RHS_3 + 
              this%BuoyancyFact = this%BulkRichardson
              call message(1,"Buoyancy term type 2 selected. Buoyancy term &
                                & calculation term uses")
              call message(2,"Bulk Richardson number:", this%BulkRichardson)
           case (3) 
              ! Parameter : Rayleigh number (avoid)  
              ! RHS_3 = RHS_3 +  
              this%BuoyancyFact = this%Ra/(this%PrandtlFluid*this%Re*this%Re) 
              call message(1,"Buoyancy term type 3 selected. Buoyancy term &
                                & calculation term uses")
              call message(2,"Rayleigh number:", this%Ra)
              call message(2,"Reynolds number:", this%Re)
           end select

           call this%sgsModel%set_BuoyancyFactor(this%BuoyancyFact)
       elseif (this%initSpinup) then
            this%BuoyancyFact = one/(this%Fr*this%Fr*this%ThetaRef)
            call this%sgsModel%set_BuoyancyFactor(this%BuoyancyFact)
       end if 

       ! STEP 22a: Set moisture
       if(this%useMoisture) then
         if(this%usescalars) then
            this%n_scalars = this%n_scalars + 1
            this%moistureIndex = this%n_scalars
         else
            if (this%useMoisture) this%usescalars = .true.
            this%n_scalars = 1
            this%moistureIndex = 1
         endif
       endif
       
       ! STEP 22b: Set other scalars
       if (this%usescalars) then
           if (allocated(this%scalars)) deallocate(this%scalars)
           allocate(this%scalars(this%n_scalars))
           do idx = 1,this%n_scalars
              if(this%useMoisture .and. idx==this%n_scalars) then
                scalar_info_dir = moisture_info_dir
              endif
              call this%scalars(idx)%init(this%gpC,this%gpE,this%spectC,this%spectE,this%sgsmodel,this%Pade6opZ,&
                           & inputfile,scalar_info_dir,this%mesh,this%u,this%v,this%w,this%wC, this%duidxjC, this%rbuffxC, &
                           & this%rbuffyC,this%rbuffzC,this%rbuffxE,this%rbuffyE,this%rbuffzE,  &
                           & this%cbuffyC,this%cbuffzC,this%cbuffyE,this%cbuffzE, this%Re, &
                           & this%isinviscid, this%useSGS, idx, this%inputdir, this%outputdir, &
                           & this%runID, useRestartFile, restartfile_TID, this%usefringe, this%usedoublefringex, &
                           & this%fringe_x, this%fringe_x1, this%fringe_x2)
           end do 
           call message(0, "SCALAR fields initialized successfully.")
       end if  
         
       ! STEP 23: Compute Rapid and Slow Pressure Split
       if (this%computeRapidSlowPressure) then
           if (this%computeDNSpressure) then
               call this%initialize_Rapid_Slow_Pressure_Split(MeanTIDX, MeanRID, MeanFilesDir)
           else
               call gracefulExit("Rapid and Slow pressure calculations require calculation of DNS pressure",13)
           end if 
       end if 
       
       ! STEP 24: Compute pressure (REDACTED) 
       !if ((this%storePressure) .or. (this%fastCalcPressure)) then
       !    call this%ComputePressure()
       !end if 

       ! STEP 25: Schedule time dumps
       this%vizDump_Schedule = vizDump_Schedule
       this%DumpThisStep = .false. 
       if (this%vizDump_Schedule == 1) then
           this%deltaT_dump = deltaT_dump
           if (useRestartFile) then
               this%t_NextDump = this%tsim - mod(this%tsim,deltaT_dump) + deltaT_dump
           else
               this%t_NextDump = this%tsim + deltaT_dump
           end if 
       end if 

       ! STEP 26: HDF5 IO
       if (ioType .ne. 0) then
           call this%initialize_hdf5_io()
           call message(0, "HDF5 IO successfully initialized.")
       end if 

       ! STEP 27: Frame angle controller 
       !! Set angle control PI yaw
       if (this%useControl) then
              allocate(this%angCont_yaw)
              call this%angCont_yaw%init(inputfile, this%spectC, this%spectE, this%gpC, this%gpE, & 
                       this%rbuffxC, this%rbuffxE, this%cbuffyC, this%cbuffyE, & 
                       this%rbuffyC, this%rbuffzC, this%restartPhi) 
       end if
       this%angleHubHeight = 1.d0       
       this%totalAngle = 0.d0
       this%wFilt = 0.d0
       this%deltaGalpha = 0.d0

       ! STEP 28: Compute the timestep
       call this%compute_deltaT()
       this%dtOld = this%dt
       this%dtRat = one 
      

       ! STEP 29: Set the buoyancy direction
       this%buoyancyDirection = buoyancyDirection
       this%useforcedStratification = useforcedStratification


       ! STEP 30: Safeguard against user invalid user inputs
       if ((this%vizDump_Schedule == 1) .and. (.not. this%useCFL)) then
           call GracefulExit("Cannot use vizDump_Schedule=1 if using fixed dt.",123)
       end if 
       if ((this%fastCalcPressure) .and. ((TimeSteppingScheme .ne. 1) .and. (TimeSteppingScheme .ne. 2))) then
           !call GracefulExit("fastCalcPressure feature is only supported with TVD RK3 or SSP RK45 time stepping.",123)
       end if

       if ((this%usescalars) .and. ((TimeSteppingScheme .ne. 1) .and. (TimeSteppingScheme .ne. 2))) then
           call GracefulExit("SCALARS are only supported with TVD RK3 or SSP RK45 time stepping.",123)
       end if 

       if ((this%fastCalcPressure) .and. (useDealiasFilterVert)) then
           call GracefulExit("fastCalcPressure feature is not supported if useDealiasFilterVert is TRUE",123) 
       end if 

       if ((.not. this%fastCalcPressure) .and. ((this%computefringePressure) .or. (this%computeDNSPressure))) then
           call gracefulExit("You need to set FASTCALCPRESSURE = .true. in & 
                    & order to use computefringepressure or computeDNSpressure", 313)
       end if
 
       if ((this%isStratified .or. this%initspinup) .and. (.not. ComputeStokesPressure )) then
          call GracefulExit("You must set ComputeStokesPressure to TRUE if &
          & there is stratification in the problem", 323)
       end if

       if (this%donot_dealias) then
           call message(0,"DONOT_DEALIAS set to TRUE. Screw you.")
       end if 
     
       if ((this%useMoisture) .and. (.not. this%isStratified)) then
           call GracefulExit("isStratified must be true for moisture to be true",123)
       end if 

       call message("IGRID initialized successfully!")
       call message("===========================================================")

   end subroutine

   
   subroutine timeAdvance(this, dtforced)
       class(igrid), intent(inout) :: this
       real(rkind), intent(in), optional :: dtforced

       select case (this%timeSteppingScheme)
       case(0)
           call this%AdamsBashforth()
       case(1)
           if(present(dtforced)) then
             call this%TVD_rk3(dtforced)
           else
             call this%TVD_rk3()
           endif
        case(2) 
           if(present(dtforced)) then
             call this%SSP_rk45(dtforced)
           else
             call this%SSP_rk45()
           endif
        case(3)
           ! Does the exact same operations as case 2
           ! written to debug case 2
           if(present(dtforced)) then
             call this%advance_SSP_RK45_all_stages(dtforced)
           else
             call this%advance_SSP_RK45_all_stages()
           end if
        end select

   end subroutine


   subroutine destroy(this)
       class(igrid), intent(inout) :: this
       integer :: idx

       if(this%useHITForcing) then
         call this%hitforce%destroy()
         deallocate(this%hitforce)
       endif 
       if (this%timeAvgFullFields) then
           call this%finalize_stats3d()
       else
       !    call this%finalize_stats()
       end if 
       nullify(this%u, this%uhat, this%v, this%vhat, this%w, this%what, this%wC)
       deallocate(this%PfieldsC, this%PfieldsE, this%SfieldsC, this%SfieldsE)
       nullify(this%u_rhs, this%v_rhs, this%w_rhs)
       deallocate(this%rhsC, this%rhsE, this%OrhsC, this%OrhsE)
       deallocate(this%duidxjC, this%duidxjChat)
       call this%spectC%destroy()
       call this%spectE%destroy()
       deallocate(this%spectC, this%spectE)
       nullify(this%nu_SGS, this%c_SGS, this%tauSGS_ij)
       if (this%useSGS) then
          call this%sgsModel%destroy()
          deallocate(this%sgsModel)
       end if

       if (allocated(this%scalars)) then
           do idx = 1,this%n_scalars
              call this%scalars(idx)%destroy()
           end do 
           deallocate(this%scalars)
       end if
   end subroutine
   
   subroutine printDivergence(this)
      class(igrid), intent(inout) :: this
      call this%padepoiss%DivergenceCheck(this%uhat, this%vhat, this%what, this%divergence)
   end subroutine 


   subroutine initialize_Rapid_Slow_Pressure_Split(this, MeanTIDX, MeanRID, MeanFilesDir)
       class(igrid), intent(inout) :: this
       integer, intent(in) :: MeanTIDX, MeanRID
       character(len=*), intent(in) :: MeanFilesDir
       real(rkind) :: maxdiv 

       if (.not. this%periodicInZ) then 
            call GracefulExit("Rapid Pressure calculation only supported for fully periodic problems",34)
       end if 

       call message(0,"Initializing the rapid/slow pressure decompositions")
       allocate(this%uM(this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3)))
       allocate(this%vM(this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3)))
       allocate(this%wM(this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3)))
       
       allocate(this%dumdx(this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3)))
       allocate(this%dvmdx(this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3)))
       allocate(this%dwmdx(this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3)))
       
       allocate(this%dumdy(this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3)))
       allocate(this%dvmdy(this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3)))
       allocate(this%dwmdy(this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3)))
       
       allocate(this%dumdz(this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3)))
       allocate(this%dvmdz(this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3)))
       allocate(this%dwmdz(this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3)))
       
       allocate(this%prapid(this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3)))
       allocate(this%pslow(this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3)))

       call readField3D(MeanRID,MeanTIDX, MeanFilesDir, "uMmn", this%uM,this%gpC)
       call readField3D(MeanRID,MeanTIDX, MeanFilesDir, "vMmn", this%vM,this%gpC)
       call readField3D(MeanRID,MeanTIDX, MeanFilesDir, "wMmn", this%wM,this%gpC)

       call this%spectC%fft(this%uM,this%cbuffyC(:,:,:,1))
       call this%spectC%mtimes_ik1_oop(this%cbuffyC(:,:,:,1),this%cbuffyC(:,:,:,2))
       call this%spectC%ifft(this%cbuffyC(:,:,:,2), this%dumdx)
       call this%spectC%mtimes_ik2_oop(this%cbuffyC(:,:,:,1),this%cbuffyC(:,:,:,2))
       call this%spectC%ifft(this%cbuffyC(:,:,:,2), this%dumdy)
       call transpose_x_to_y(this%uM, this%rbuffyC(:,:,:,1), this%gpC)
       call transpose_y_to_z(this%rbuffyC(:,:,:,1), this%rbuffzC(:,:,:,1), this%gpC)
       call this%Pade6opZ%ddz_C2C(this%rbuffzC(:,:,:,1),this%rbuffzC(:,:,:,2), uBC_bottom, uBC_top)
       call transpose_z_to_y(this%rbuffzC(:,:,:,2),this%rbuffyC(:,:,:,1), this%gpC)
       call transpose_y_to_x(this%rbuffyC(:,:,:,1),this%dumdz, this%gpC)

       call this%spectC%fft(this%vM,this%cbuffyC(:,:,:,1))
       call this%spectC%mtimes_ik1_oop(this%cbuffyC(:,:,:,1),this%cbuffyC(:,:,:,2))
       call this%spectC%ifft(this%cbuffyC(:,:,:,2), this%dvmdx)
       call this%spectC%mtimes_ik2_oop(this%cbuffyC(:,:,:,1),this%cbuffyC(:,:,:,2))
       call this%spectC%ifft(this%cbuffyC(:,:,:,2), this%dvmdy)
       call transpose_x_to_y(this%vM, this%rbuffyC(:,:,:,1), this%gpC)
       call transpose_y_to_z(this%rbuffyC(:,:,:,1), this%rbuffzC(:,:,:,1), this%gpC)
       call this%Pade6opZ%ddz_C2C(this%rbuffzC(:,:,:,1),this%rbuffzC(:,:,:,2), vBC_bottom, vBC_top)
       call transpose_z_to_y(this%rbuffzC(:,:,:,2),this%rbuffyC(:,:,:,1), this%gpC)
       call transpose_y_to_x(this%rbuffyC(:,:,:,1),this%dvmdz, this%gpC)

       call this%spectC%fft(this%wM,this%cbuffyC(:,:,:,1))
       call this%spectC%mtimes_ik1_oop(this%cbuffyC(:,:,:,1),this%cbuffyC(:,:,:,2))
       call this%spectC%ifft(this%cbuffyC(:,:,:,2), this%dwmdx)
       call this%spectC%mtimes_ik2_oop(this%cbuffyC(:,:,:,1),this%cbuffyC(:,:,:,2))
       call this%spectC%ifft(this%cbuffyC(:,:,:,2), this%dwmdy)
       call transpose_x_to_y(this%wM, this%rbuffyC(:,:,:,1), this%gpC)
       call transpose_y_to_z(this%rbuffyC(:,:,:,1), this%rbuffzC(:,:,:,1), this%gpC)
       call this%Pade6opZ%ddz_C2C(this%rbuffzC(:,:,:,1),this%rbuffzC(:,:,:,2), wBC_bottom, wBC_top)
       call transpose_z_to_y(this%rbuffzC(:,:,:,2),this%rbuffyC(:,:,:,1), this%gpC)
       call transpose_y_to_x(this%rbuffyC(:,:,:,1),this%dwmdz, this%gpC)
        
       this%rbuffxC(:,:,:,1) = abs(this%dumdx + this%dvmdy + this%dwmdz)
       maxdiv = p_maxval(maxval(this%rbuffxC(:,:,:,1)))

       call this%poiss_periodic%init(this%dx, this%dy, this%dz, this%gpC, 1, useExhaustiveFFT=.true.)
       
       call message(1,"Read in the mean fields and gradients needed for pressure decompositions")
       call message(1,"Maximum divergence in mean fields:", maxdiv)
   end subroutine 
   
   
   
   subroutine get_boundary_conditions_stencil()

        wBC_bottom     = -1; wBC_top     = -1;  
        WdWdzBC_bottom = -1; WdWdzBC_top = -1;
        WWBC_bottom    = +1; WWBC_top    = +1;
        dWdzBC_bottom  =  0; dWdzBC_top  =  0;

        !! Bottom wall 
        call message(0,"Bottom Wall Boundary Condition is:")
        select case (botWall)
        case(1)
           call message(1,"No-Slip Wall")
           ! NOTE: no-slip wall requires both w = 0 and dwdz = 0. Therefore, w
           ! is an even extension, which also satisfies w = 0.
           !uBC_bottom      = -1; vBC_bottom      = -1;
           uBC_bottom      =  0; vBC_bottom      =  0;
           dUdzBC_bottom   =  0; dVdzBC_bottom   =  0;
           WdUdzBC_bottom  =  0; WdVdzBC_bottom  =  0;
           !UWBC_bottom     = +1; VWBC_bottom     = +1;
           UWBC_bottom     =  0; VWBC_bottom     =  0;
           wBC_bottom      = +1; WdWdzBC_bottom  = -1; 
           WWBC_bottom     = +1; dwdzBC_bottom   = -1;
        case(2) 
           call message(1,"Slip Wall")
           uBC_bottom      = +1; vBC_bottom      = +1;
           dUdzBC_bottom   = -1; dVdzBC_bottom   = -1;
           WdUdzBC_bottom  = +1; WdVdzBC_bottom  = +1;
           UWBC_bottom     = -1; VWBC_bottom     = -1;   
        case(3) 
           call message(1,"Wall Model")
           uBC_bottom      =  0; vBC_bottom      =  0;
           dUdzBC_bottom   =  0; dVdzBC_bottom   =  0;
           WdUdzBC_bottom  = -1; WdVdzBC_bottom  = -1;
           UWBC_bottom     = -1; VWBC_bottom     = -1;   
        case default
           call gracefulExit("Invalid choice for BOTTOM WALL BCs",423)
        end select
        
        ! Top wall 
        call message(0,"Top Wall Boundary Condition is:")
        select case (TopWall)
        case(1)
           call message(1,"No-Slip Wall")
           ! NOTE: no-slip wall requires both w = 0 and dwdz = 0. Therefore, w
           ! is an even extension, which also satisfies w = 0.
           ! uBC_top      = -1; vBC_top      = -1;
           uBC_top      =  0; vBC_top      =  0;
           dUdzBC_top   =  0; dVdzBC_top   =  0;
           WdUdzBC_top  =  0; WdVdzBC_top  =  0;
           !UWBC_top     = +1; VWBC_top     = +1;   
           UWBC_top     =  0; VWBC_top     =  0;   
           wBC_top      = +1; WdWdzBC_top  = -1; 
           WWBC_top     = +1; dwdzBC_top   = -1;
        case(2) 
           call message(1,"Slip Wall")
           uBC_top      = +1; vBC_top      = +1;
           dUdzBC_top   = -1; dVdzBC_top   = -1;
           WdUdzBC_top  = +1; WdVdzBC_top  = +1;
           UWBC_top     = -1; VWBC_top     = -1;   
        case(3) 
           call message(1,"Wall Model")
           uBC_top      =  0; vBC_top      =  0;
           dUdzBC_top   =  0; dVdzBC_top   =  0;
           WdUdzBC_top  = -1; WdVdzBC_top  = -1;
           UWBC_top     = -1; VWBC_top     = -1;   
        case default
           call gracefulExit("Invalid choice for TOP WALL BCs",13)
        end select
       
        select case (topBC_Temp)
        case(0) ! Dirichlet (default)
           TBC_top = 0; dTdzBC_top = 0; WTBC_top = -1;
           WdTdzBC_top = 0;
        case(1)
           TBC_top = 1; dTdzBC_top = -1; WTBC_top = -1;
           WdTdzBC_top = 1;
         case (2) ! Inhomogeneous Neumann BC for temperature at the top
            TBC_top = 0; dTdzBC_top = 0; WTBC_top = -1;
            WdTdzBC_top = 0;
        case (3)
           TBC_top = 0; dTdzBC_top = 0; WTBC_top = -1;
           WdTdzBC_top = 0;
        end select 
        select case (botBC_Temp)
        case (0) ! Dirichlet BC for temperature at the bottom
           TBC_bottom = 0; dTdzBC_bottom = 0; WTBC_bottom = -1; 
           WdTdzBC_bottom = 0;      
        case(1)  ! Homogenenous Neumann BC at the bottom
           TBC_bottom = 1; dTdzBC_bottom = -1; WTBC_bottom = -1;
           WdTdzBC_bottom = 1
         case (2) ! Inhomogeneous Neumann BC for temperature at the bottom
            TBC_bottom = 0; dTdzBC_bottom = 0; WTBC_bottom = -1;
            WdTdzBC_bottom = 0;
        case (3) 
           TBC_bottom = 0; dTdzBC_bottom = 0; WTBC_bottom = -1; 
           WdTdzBC_bottom = 0;      
        end select

   end subroutine


   subroutine initialize_scalar_for_initspinup(this,useRestartFile, inputfile, tid, rid)
       use random,             only: gaussian_random
       use decomp_2d_io
       class(igrid), intent(inout) :: this 
       real(rkind), dimension(:,:,:), allocatable :: randArr
       real(rkind)  :: ScalarRef = 290.d0, sigma_purt = 0.1d0 
       real(rkind)  :: Froude_Scalar = 0.08d0
       real(rkind)  :: zcutoff = 0.25d0, Tstop_initSpinUp = 10.d0
       logical, intent(in) :: useRestartFile
       character(len=clen), intent(in) :: inputfile
       character(len=clen) :: fname, tempname
       integer, intent(in) :: tid, rid
       integer :: seed = 2123122

       namelist /INIT_SPINUP/ Tstop_InitSpinUp, Zcutoff, ScalarRef, Sigma_purt, Froude_Scalar, seed

       open(unit=11, file=trim(inputfile), form='FORMATTED', iostat=ierr)
       read(unit=11, NML=INIT_SPINUP)
       close(11)

       this%Tstop_InitSpinUp = Tstop_InitSpinUp 
       this%Fr = Froude_Scalar
       this%ThetaRef = ScalarRef
       if (useRestartFile) then
           write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",rid, "_T.",tid
           fname = this%InputDir(:len_trim(this%InputDir))//"/"//trim(tempname)
           call decomp_2d_read_one(1,this%T,fname, this%gpC)
           call message(0,"Read the spinup scalar field from the restart file")
       else
           allocate(randArr(size(this%T,1),size(this%T,2),size(this%T,3)))
           this%T = ScalarRef
           call gaussian_random(randArr,zero,one,seed + 10*nrank)
           randArr = sigma_purt*randArr
           ! Set random numbers in z > 0.25 to zero 
           where (this%mesh(:,:,:,3) > Zcutoff)
               randArr = zero
           end where
           
           ! Add Temperature purturbations 
           this%T = this%T + randArr
           call message(0,"Initialized the spinup scalar field") 
           deallocate(randArr)
       end if

       TBC_bottom = +1; dTdzBC_bottom = -1; WTBC_bottom = -1;
       WdTdzBC_bottom = +1
       TBC_top = +1; dTdzBC_top = -1; WTBC_top = -1;
       WdTdzBC_top = +1

       if (this%tsim > Tstop_InitSpinUp) this%initspinup = .false. 
   end subroutine 

   subroutine getMaxOddballModes(this,&
                                 uXoddMaxI, uYoddMaxI, &
                                 uXoddMaxR, uYoddMaxR, &
                                 vXoddMaxI, vYoddMaxI, &
                                 vXoddMaxR, vYoddMaxR, &
                                 wXoddMaxI, wYoddMaxI, &
                                 wXoddMaxR, wYoddMaxR, dataLoc)
       use fortran_assert, only: assert

       class(igrid), intent(inout) :: this
       real(rkind), intent(out) :: uXoddMaxI, uYoddMaxI, &
                                   uXoddMaxR, uYoddMaxR, & 
                                   vXoddMaxI, vYoddMaxI, & 
                                   vXoddMaxR, vYoddMaxR, & 
                                   wXoddMaxI, wYoddMaxI, & 
                                   wXoddMaxR, wYoddMaxR 
       character(len=1), intent(in) :: dataLoc
       integer :: nx, ny, nz
       complex(rkind), dimension(:,:,:,:), allocatable :: cbuffxC, cbuffxE

       nx = this%gpC%xsz(1)
       ny = this%gpC%ysz(2)
       nz = this%gpC%zsz(3)
       
       select case (dataLoc)
       case ('C')
         allocate(cbuffxC(this%sp_gpC%xsz(1), this%sp_gpC%xsz(2), &
           this%sp_gpC%xsz(3), 3))
         uYoddMaxI = p_maxval(maxval(aimag(this%uhat (:,ny/2+1,:))))
         vYoddMaxI = p_maxval(maxval(aimag(this%vhat (:,ny/2+1,:))))
         wYoddMaxI = p_maxval(maxval(aimag(this%whatC(:,ny/2+1,:))))

         uYoddMaxR = p_maxval(maxval(real(this%uhat (:,ny/2+1,:),rkind)))
         vYoddMaxR = p_maxval(maxval(real(this%vhat (:,ny/2+1,:),rkind)))
         wYoddMaxR = p_maxval(maxval(real(this%whatC(:,ny/2+1,:),rkind)))
         
         call transpose_y_to_x(this%uhat, cbuffxC(:,:,:,1),this%sp_gpC)
         call transpose_y_to_x(this%vhat, cbuffxC(:,:,:,2),this%sp_gpC)
         call transpose_y_to_x(this%whatC,cbuffxC(:,:,:,3),this%sp_gpC)
         
         uXoddMaxI = p_maxval(maxval(aimag(cbuffxC(nx/2+1,:,:,1))))
         vXoddMaxI = p_maxval(maxval(aimag(cbuffxC(nx/2+1,:,:,2))))
         wXoddMaxI = p_maxval(maxval(aimag(cbuffxC(nx/2+1,:,:,3))))
         
         uXoddMaxR = p_maxval(maxval(real(cbuffxC(nx/2+1,:,:,1),rkind)))
         vXoddMaxR = p_maxval(maxval(real(cbuffxC(nx/2+1,:,:,2),rkind)))
         wXoddMaxR = p_maxval(maxval(real(cbuffxC(nx/2+1,:,:,3),rkind)))
         
         deallocate(cbuffxC)

       case ('E')
         allocate(cbuffxE(this%sp_gpE%xsz(1), this%sp_gpE%xsz(2), &
           this%sp_gpE%xsz(3), 3))
         uYoddMaxI = p_maxval(maxval(aimag(this%uEhat(:,ny/2+1,:))))
         vYoddMaxI = p_maxval(maxval(aimag(this%vEhat(:,ny/2+1,:))))
         wYoddMaxI = p_maxval(maxval(aimag(this%what (:,ny/2+1,:))))

         uYoddMaxR = p_maxval(maxval(real(this%uEhat(:,ny/2+1,:),rkind)))
         vYoddMaxR = p_maxval(maxval(real(this%vEhat(:,ny/2+1,:),rkind)))
         wYoddMaxR = p_maxval(maxval(real(this%what (:,ny/2+1,:),rkind)))

         call transpose_y_to_x(this%uEhat,cbuffxE(:,:,:,1),this%sp_gpE)
         call transpose_y_to_x(this%vEhat,cbuffxE(:,:,:,2),this%sp_gpE)
         call transpose_y_to_x(this%what, cbuffxE(:,:,:,3),this%sp_gpE)
         
         uXoddMaxI = p_maxval(maxval(aimag(cbuffxE(nx/2+1,:,:,1))))
         vXoddMaxI = p_maxval(maxval(aimag(cbuffxE(nx/2+1,:,:,2))))
         wXoddMaxI = p_maxval(maxval(aimag(cbuffxE(nx/2+1,:,:,3))))

         uXoddMaxR = p_maxval(maxval(real(cbuffxE(nx/2+1,:,:,1),rkind)))
         vXoddMaxR = p_maxval(maxval(real(cbuffxE(nx/2+1,:,:,2),rkind)))
         wXoddMaxR = p_maxval(maxval(real(cbuffxE(nx/2+1,:,:,3),rkind)))
         
         deallocate(cbuffxE)
       end select

   end subroutine

end module 
