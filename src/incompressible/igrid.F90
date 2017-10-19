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
    use reductions, only: p_maxval, p_sum
    use timer, only: tic, toc
    use PadePoissonMod, only: Padepoisson 
    use sgsmod_igrid, only: sgs_igrid
    use wallmodelMod, only: wallmodel
    use numerics
    !use cd06staggstuff, only: cd06stagg
    use cf90stuff, only: cf90
    use TurbineMod, only: TurbineArray 
    use kspreprocessing, only: ksprep  
    use PadeDerOps, only: Pade6Stagg
    use Fringemethod, only: fringe
    use forcingmod,   only: HIT_shell_forcing

    implicit none

    private
    public :: igrid 

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
        real(rkind) :: dt, tstop, CFL, dx, dy, dz, tsim
        character(len=clen) ::  outputdir
        real(rkind), dimension(:,:,:,:), allocatable :: mesh
        character(len=clen) :: filter_x          ! What filter to use in X: "cf90", "gaussian", "lstsq", "spectral"
        character(len=clen) :: filter_y          ! What filter to use in X: "cf90", "gaussian", "lstsq", "spectral" 
        character(len=clen) :: filter_z          ! What filter to use in X: "cf90", "gaussian", "lstsq", "spectral" 
        integer      :: step, nsteps = 999999

        type(decomp_info), allocatable :: gpC, gpE
        type(decomp_info), pointer :: Sp_gpC, Sp_gpE
        type(spectral), allocatable :: spectE, spectC
        type(staggOps), allocatable :: Ops, OpsPP
        type(sgs_igrid), allocatable :: sgsmodel
        type(wallmodel), allocatable :: moengWall

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
        real(rkind), dimension(:,:,:), allocatable :: dTdxC, dTdyC, dTdzE, dTdzC

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
    
        integer :: BuoyancyTermType = 0 
        real(rkind) :: BuoyancyFact = 0.d0

        logical :: periodicInZ  = .false. 
        logical :: newTimeStep = .true., computeVorticity = .false.  
        integer :: timeSteppingScheme = 0 
        integer :: runID, t_start_planeDump, t_stop_planeDump, t_planeDump, t_DivergenceCheck
        integer :: t_start_pointProbe, t_stop_pointProbe, t_pointProbe
        logical :: useCoriolis = .true. , isStratified = .false., useSponge = .false. 
        logical :: useExtraForcing = .false., useGeostrophicForcing = .false., isInviscid = .false.  
        logical :: useSGS = .false. 
        logical :: UseDealiasFilterVert = .false.
        logical :: useDynamicProcedure 
        logical :: useCFL = .false.  
        logical :: dumpPlanes = .false., useWindTurbines = .false. 

        complex(rkind), dimension(:,:,:), allocatable :: dPf_dxhat

        real(rkind) :: max_nuSGS, invObLength, Tsurf0, Tsurf, dTsurf_dt, ThetaRef

        real(rkind) :: dtOld, dtRat, Tmn, wTh_surf
        real(rkind), dimension(:,:,:), allocatable :: filteredSpeedSq
        integer :: wallMType, botBC_Temp 

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
        logical :: computeForcingTerm = .false. 

        ! Pointers linked to SGS stuff
        real(rkind), dimension(:,:,:,:), pointer :: tauSGS_ij
        real(rkind), dimension(:,:,:)  , pointer :: nu_SGS, tau13, tau23
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
        real(rkind) :: KSFilFact 
        real(rkind), dimension(:,:,:), allocatable :: KS_probe_data


        ! Pressure Solver
        logical :: StorePressure = .false., fastCalcPressure = .true. 
        integer :: P_dumpFreq = 10, P_compFreq = 10
        logical :: AlreadyHaveRHS = .false.
        real(rkind), dimension(:,:,:), allocatable :: pressure  

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

        ! HIT Forcing
        logical :: useHITForcing = .false.
        type(HIT_shell_forcing), allocatable :: hitforce

        contains
            procedure          :: init
            procedure          :: destroy
            procedure          :: printDivergence 
            procedure          :: getMaxKE
            procedure          :: getMeanKE
            procedure          :: timeAdvance
            procedure          :: start_io
            procedure          :: finalize_io
            procedure          :: get_dt
            procedure          :: interpolate_cellField_to_edgeField 
            !procedure, private :: init_stats
            procedure, private :: debug
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
            procedure, private :: dumpRestartFile
            procedure, private :: readRestartFile
            procedure, private :: compute_z_mean 
            procedure, private :: compute_deltaT
            procedure, private :: dump_pointProbes
            procedure, private :: dump_stats3D
            procedure, private :: compute_stats3D 
            procedure, private :: dealiasFields
            procedure, private :: ApplyCompactFilter 
            procedure, private :: addNonLinearTerm_skewSymm
            procedure, private :: populate_rhs
            procedure, private :: project_and_prep
            procedure, private :: wrapup_timestep
            procedure, private :: reset_pointers
            procedure, private :: compute_vorticity
            procedure, private :: finalize_stats3D
            procedure, private :: dump_planes
            procedure          :: dumpFullField 
            procedure, private :: dumpVisualizationInfo
            procedure, private :: DeletePrevStats3DFiles
            procedure, private :: Delete_file_if_present
            procedure, private :: updateProbes 
            procedure, private :: dumpProbes
            procedure, private :: correctPressureRotationalForm
            procedure, private :: initialize_scalar_for_InitSpinUp
    end type

contains 

    subroutine init(this,inputfile, initialize2decomp)
        class(igrid), intent(inout), target :: this        
        character(len=clen), intent(in) :: inputfile 
        character(len=clen) :: outputdir, inputdir, turbInfoDir, ksOutputDir, controlDir = "null"
        integer :: nx, ny, nz, prow = 0, pcol = 0, ioUnit, nsteps = 999999
        integer :: tid_StatsDump =10000, tid_compStats = 10000,  WallMType = 0, t_planeDump = 1000
        integer :: t_pointProbe = 10000, t_start_pointProbe = 10000, t_stop_pointProbe = 1
        integer :: runID = 0,  t_dataDump = 99999, t_restartDump = 99999,t_stop_planeDump = 1,t_dumpKSprep = 10 
        integer :: restartFile_TID = 1, ioType = 0, restartFile_RID =1, t_start_planeDump = 1
        real(rkind) :: dt=-one,tstop=one,CFL =-one,tSimStartStats=100.d0,dpfdy=zero,dPfdz=zero,ztop
        real(rkind) :: Pr = 0.7_rkind, Re = 8000._rkind, Ro = 1000._rkind,dpFdx = zero, G_alpha = 0.d0, PrandtlFluid = 1.d0
        real(rkind) :: SpongeTscale = 50._rkind, zstSponge = 0.8_rkind, Fr = 1000.d0, G_geostrophic = 1.d0
        logical ::useRestartFile=.false.,isInviscid=.false.,useCoriolis = .true., PreProcessForKS = .false.  
        logical ::isStratified=.false.,dumpPlanes = .false.,useExtraForcing = .false.
        logical ::useSGS = .false.,useSpongeLayer=.false.,useWindTurbines = .false., useTopAndBottomSymmetricSponge = .false. 
        logical :: useGeostrophicForcing = .false., PeriodicInZ = .false. 
        real(rkind), dimension(:,:,:), pointer :: zinZ, zinY, zEinY, zEinZ
        integer :: AdvectionTerm = 1, NumericalSchemeVert = 0, t_DivergenceCheck = 10, ksRunID = 10
        integer :: timeSteppingScheme = 0, num_turbines = 0, P_dumpFreq = 10, P_compFreq = 10, BuoyancyTermType = 1
        logical :: normStatsByUstar=.false., ComputeStokesPressure = .true., UseDealiasFilterVert = .false.
        real(rkind) :: Lz = 1.d0, latitude = 90._rkind, KSFilFact = 4.d0, dealiasFact = 2.d0/3.d0, frameAngle = 0.d0, BulkRichardson = 0.d0
        logical :: ADM = .false., storePressure = .false., useSystemInteractions = .true., useFringe = .false., useHITForcing = .false.
        integer :: tSystemInteractions = 100, ierr, KSinitType = 0, nKSvertFilt = 1, ADM_Type = 1
        logical :: computeSpectra = .false., timeAvgFullFields = .false., fastCalcPressure = .true., usedoublefringex = .false.  
        logical :: assume_fplane = .true., periodicbcs(3), useProbes = .false., KSdoZfilter = .true., computeVorticity = .false.  
        real(rkind), dimension(:,:), allocatable :: probe_locs
        real(rkind), dimension(:), allocatable :: temp
        integer :: ii, idx, temploc(1)
        logical, intent(in), optional :: initialize2decomp
        logical :: reset2decomp, InitSpinUp = .false., useExhaustiveFFT = .true.  

        namelist /INPUT/ nx, ny, nz, tstop, dt, CFL, nsteps, inputdir, outputdir, prow, pcol, &
                         useRestartFile, restartFile_TID, restartFile_RID 
        namelist /IO/ t_restartDump, t_dataDump, ioType, dumpPlanes, runID, useProbes, &
                        t_planeDump, t_stop_planeDump, t_start_planeDump, t_start_pointProbe, t_stop_pointProbe, t_pointProbe
        namelist /STATS/tid_StatsDump,tid_compStats,tSimStartStats,normStatsByUstar,computeSpectra,timeAvgFullFields, computeVorticity
        namelist /PHYSICS/isInviscid,useCoriolis,useExtraForcing,isStratified,Re,Ro,Pr,Fr, useSGS, PrandtlFluid, BulkRichardson, BuoyancyTermType,&
                          useGeostrophicForcing, G_geostrophic, G_alpha, dpFdx, dpFdy, dpFdz, assume_fplane, latitude, useHITForcing, frameAngle
        namelist /BCs/ PeriodicInZ, topWall, botWall, useSpongeLayer, zstSponge, SpongeTScale, botBC_Temp, topBC_Temp, useTopAndBottomSymmetricSponge, useFringe, usedoublefringex
        namelist /WINDTURBINES/ useWindTurbines, num_turbines, ADM, turbInfoDir, ADM_Type  
        namelist /NUMERICS/ AdvectionTerm, ComputeStokesPressure, NumericalSchemeVert, &
                            UseDealiasFilterVert, t_DivergenceCheck, TimeSteppingScheme, InitSpinUp, &
                                 useExhaustiveFFT, dealiasFact 
        namelist /KSPREPROCESS/ PreprocessForKS, KSoutputDir, KSRunID, t_dumpKSprep, KSinitType, KSFilFact, &
                                 KSdoZfilter, nKSvertFilt
        namelist /PRESSURE_CALC/ fastCalcPressure, storePressure, P_dumpFreq, P_compFreq            
        namelist /OS_INTERACTIONS/ useSystemInteractions, tSystemInteractions, controlDir

        ! STEP 1: READ INPUT 
        ioUnit = 11
        open(unit=ioUnit, file=trim(inputfile), form='FORMATTED', iostat=ierr)
        read(unit=ioUnit, NML=INPUT)
        read(unit=ioUnit, NML=NUMERICS)
        read(unit=ioUnit, NML=IO)
        read(unit=ioUnit, NML=STATS)
        read(unit=ioUnit, NML=OS_INTERACTIONS)
        read(unit=ioUnit, NML=PHYSICS)
        read(unit=ioUnit, NML=PRESSURE_CALC)
        read(unit=ioUnit, NML=BCs)
        read(unit=ioUnit, NML=WINDTURBINES)
        read(unit=ioUnit, NML=KSPREPROCESS)
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
        this%KSinitType = KSinitType; this%KSFilFact = KSFilFact; this%useFringe = useFringe
        this%nsteps = nsteps; this%PeriodicinZ = periodicInZ; this%usedoublefringex = usedoublefringex 
        this%useHITForcing = useHITForcing; this%BuoyancyTermType = BuoyancyTermType 
        this%frameAngle = frameAngle; this%computeVorticity = computeVorticity

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
        Lz = p_maxval(this%mesh(:,:,:,3)) + this%dz/2.d0
        call message(0,"Mesh generated:")
        call message(1,"dx:", this%dx)
        call message(1,"dy:", this%dy)
        call message(1,"dz:", this%dz)
        call message(1,"Lz:", Lz)


        ! STEP 4: ALLOCATE/INITIALIZE THE SPECTRAL DERIVED TYPES
        allocate(this%spectC)
        call this%spectC%init("x", nx, ny, nz  , this%dx, this%dy, this%dz, &
                "four", this%filter_x, 2 , fixOddball=.false., exhaustiveFFT=useExhaustiveFFT, init_periodicInZ=periodicinZ, dealiasF=dealiasfact)
        allocate(this%spectE)
        call this%spectE%init("x", nx, ny, nz+1, this%dx, this%dy, this%dz, &
                "four", this%filter_x, 2 , fixOddball=.false., exhaustiveFFT=useExhaustiveFFT, init_periodicInZ=.false., dealiasF=dealiasfact)
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
        allocate(this%cbuffyC(this%sp_gpC%ysz(1),this%sp_gpC%ysz(2),this%sp_gpC%ysz(3),2))
        allocate(this%cbuffyE(this%sp_gpE%ysz(1),this%sp_gpE%ysz(2),this%sp_gpE%ysz(3),2))
        
        allocate(this%cbuffzC(this%sp_gpC%zsz(1),this%sp_gpC%zsz(2),this%sp_gpC%zsz(3),3))
        allocate(this%cbuffzE(this%sp_gpE%zsz(1),this%sp_gpE%zsz(2),this%sp_gpE%zsz(3),2))

        allocate(this%rbuffxC(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3),2))
        allocate(this%rbuffyC(this%gpC%ysz(1),this%gpC%ysz(2),this%gpC%ysz(3),2))
        allocate(this%rbuffzC(this%gpC%zsz(1),this%gpC%zsz(2),this%gpC%zsz(3),4))

        allocate(this%rbuffxE(this%gpE%xsz(1),this%gpE%xsz(2),this%gpE%xsz(3),2))
        allocate(this%rbuffyE(this%gpE%ysz(1),this%gpE%ysz(2),this%gpE%ysz(3),2))
        allocate(this%rbuffzE(this%gpE%zsz(1),this%gpE%zsz(2),this%gpE%zsz(3),4))
        allocate(this%filteredSpeedSq(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)))
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
            else
                call GraceFulExit("Only Dirichlet and homog. Neumann BCs supported for Temperature at &
                    & this time. Set botBC_Temp = 0 or 1",341)        
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

            allocate(this%SGSmodel)
            call this%sgsModel%init(this%gpC, this%gpE, this%spectC, this%spectE, this%dx, this%dy, this%dz, inputfile, &
                                    this%rbuffxE(1,1,:,1), this%mesh(1,1,:,3), this%fBody_x, this%fBody_y, this%fBody_z, &
                                    this%storeFbody,this%Pade6opZ, this%cbuffyC, this%cbuffzC, this%cbuffyE, this%cbuffzE, &
                                    this%rbuffxC, this%rbuffyC, this%rbuffzC, this%rbuffyE, this%rbuffzE, this%Tsurf, &
                                    this%ThetaRef, this%Fr, this%Re, Pr, this%isInviscid, this%isStratified, this%botBC_Temp, &
                                    this%initSpinUp)
            call this%sgsModel%link_pointers(this%nu_SGS, this%tauSGS_ij, this%tau13, this%tau23, this%q1, this%q2, this%q3)
            call message(0,"SGS model initialized successfully")
        end if 
        this%max_nuSGS = zero


        ! STEP 12: Set Sponge Layer
        if (this%useSponge) then
            allocate(this%RdampC(this%sp_gpC%ysz(1), this%sp_gpC%ysz(2), this%sp_gpC%ysz(3)))
            allocate(this%RdampE(this%sp_gpE%ysz(1), this%sp_gpE%ysz(2), this%sp_gpE%ysz(3)))
            zinY => this%rbuffyC(:,:,:,1); zinZ => this%rbuffzC(:,:,:,1)
            zEinZ => this%rbuffzE(:,:,:,1); zEinY => this%rbuffyE(:,:,:,1)
            call transpose_x_to_y(this%mesh(:,:,:,3),zinY,this%gpC)
            call transpose_y_to_z(zinY,zinZ,this%gpC)
            call this%OpsPP%InterpZ_Cell2Edge(zinZ,zEinZ,zero,zero)
            zEinZ(:,:,this%nz+1) = zEinZ(:,:,this%nz) + this%dz
            ztop = zEinZ(1,1,this%nz+1); zstSponge = zstSponge*ztop 
            call transpose_z_to_y(zEinZ,zEinY,this%gpE)
            this%RdampC = (one/SpongeTscale) * (one - cos(pi*(zinY - zstSponge) /(zTop - zstSponge)))/two
            this%RdampE = (one/SpongeTscale) * (one - cos(pi*(zEinY - zstSponge)/(zTop - zstSponge)))/two
            if (useTopAndBottomSymmetricSponge) then
               where (abs(zEinY) < zstSponge) 
                   this%RdampE = zero
               end where
               where (abs(zinY) < zstSponge) 
                   this%RdampC = zero
               end where
               if ((zEinZ(1,1,1) + zEinZ(1,1,this%nz+1))<1.d-13) then
                  call message(0,"WARNING: Computed domain is not symmetric &
                     & about z=0. You shouldn't use the symmetric sponge")
                  call MPI_BARRIER(mpi_comm_world, ierr)
                  call GracefulExit("Failed at sponge initialization",134)
               end if   
            else
               where (zEinY < zstSponge) 
                   this%RdampE = zero
               end where
               where (zinY < zstSponge) 
                   this%RdampC = zero
               end where
            end if 
            

            call this%spectC%alloc_r2c_out(this%uBase)
            call this%spectC%alloc_r2c_out(this%vBase)
            call this%spectC%alloc_r2c_out(this%TBase)
            this%rbuffxC(:,:,:,1) = this%u 
            call this%spectC%fft(this%rbuffxC(:,:,:,1),this%uBase)
            this%rbuffxC(:,:,:,1) = this%v
            call this%spectC%fft(this%rbuffxC(:,:,:,1),this%vBase)
            this%rbuffxC(:,:,:,1) = this%T
            call this%spectC%fft(this%rbuffxC(:,:,:,1),this%TBase)
            call message(0,"Sponge Layer initialized successfully")
            call message(1,"Sponge Layer active above z = ",zstSponge)
        end if 

        if (this%useWindTurbines) then
            allocate(this%WindTurbineArr)
            call this%WindTurbineArr%init(inputFile, this%gpC, this%gpE, this%spectC, this%spectE, this%cbuffyC, this%cbuffyE, this%cbuffzC, this%cbuffzE, this%mesh, this%dx, this%dy, this%dz)
        end if
        ! STEP 12: Set visualization planes for io
        call set_planes_io(this%xplanes, this%yplanes, this%zplanes)

        ! STEP 13: Compute the timestep
        call this%compute_deltaT()
        this%dtOld = this%dt
        this%dtRat = one 


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
                allocate(this%probe_data(1:6,1:this%nprobes,0:this%probeTimeLimit-1))
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
                call this%LES2KS%init(this%spectC, this%gpC, this%dx, this%dy, this%outputdir, this%RunID, this%probes, this%KSFilFact, KSdoZfilter, nKSvertFilt)
                call this%LES2KS%link_pointers(this%uFil4KS, this%vFil4KS, this%wFil4KS)
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
            call message(0, "KS Preprocessor initializaed successfully.")
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
        call mpi_barrier(mpi_comm_world, ierr)
        if (this%usedoublefringex) then
            allocate(this%fringe_x1, this%fringe_x2)
            call this%fringe_x1%init(inputfile, this%dx, this%mesh(:,1,1,1), this%dy, this%mesh(1,:,1,2), &
                                        this%spectC, this%spectE, this%gpC, this%gpE, &
                                        this%rbuffxC, this%rbuffxE, this%cbuffyC, this%cbuffyE, fringeID=1)   
            call this%fringe_x2%init(inputfile, this%dx, this%mesh(:,1,1,1), this%dy, this%mesh(1,:,1,2), &
                                        this%spectC, this%spectE, this%gpC, this%gpE, &
                                        this%rbuffxC, this%rbuffxE, this%cbuffyC, this%cbuffyE, fringeID=2)   

        else
            if (this%useFringe) then
                allocate(this%fringe_x)
                call this%fringe_x%init(inputfile, this%dx, this%mesh(:,1,1,1), this%dy, this%mesh(1,:,1,2), &
                                        this%spectC, this%spectE, this%gpC, this%gpE, &
                                        this%rbuffxC, this%rbuffxE, this%cbuffyC, this%cbuffyE)   
            end if
        end if 
        
        ! STEP 18: Set HIT Forcing
        if (this%useHITForcing) then
            allocate(this%hitforce)
            call this%hitforce%init(inputfile, this%sp_gpC, this%sp_gpE, this%spectC, this%cbuffyE(:,:,:,1), &
                           this%cbuffyC(:,:,:,1), this%cbuffzE(:,:,:,1), this%cbuffzC, this%step)
        end if
        
        ! STEP 19: Set up storage for Pressure
        if ((this%storePressure) .or. (this%fastCalcPressure)) then
            allocate(this%Pressure(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)))
            call this%ComputePressure()
            call message(1, "Done allocating storage for pressure")
        end if 

        ! STEP 20: Update the probes
        if (this%useProbes) call this%updateProbes()


        ! STEP 21: Buoyancy term type
        if (this%isStratified) then
            select case (this%BuoyancyTermType)
            case(1)
               this%BuoyancyFact = one/(this%Fr*this%Fr*this%ThetaRef)
               call message(1,"Buoyancy term type 1 selected. Buoyancy term &
                                 & calculation term uses")
               call message(2,"Froude number:", this%Fr)
               call message(2,"Reference temperature:", this%thetaRef)
            case(2)
               this%BuoyancyFact = this%BulkRichardson
               call message(1,"Buoyancy term type 2 selected. Buoyancy term &
                                 & calculation term uses")
               call message(2,"Bulk Richardson number:", this%BulkRichardson)
            end select
        end if 


        ! STEP 22: Safeguard against user invalid user inputs
        if ((this%fastCalcPressure) .and. ((TimeSteppingScheme .ne. 1) .and. (TimeSteppingScheme .ne. 2))) then
            call GracefulExit("fastCalcPressure feature is only supported with TVD RK3 or SSP RK45 time stepping.",123)
        end if 
        if ((this%fastCalcPressure) .and. (useDealiasFilterVert)) then
            call GracefulExit("fastCalcPressure feature is not supported if useDealiasFilterVert is TRUE",123) 
        end if 

        if ((this%isStratified .or. this%initspinup) .and. (.not. ComputeStokesPressure )) then
            call GracefulExit("You must set ComputeStokesPressure to TRUE if &
            & there is stratification in the problem", 323)
        end if 

      
        call message("IGRID initialized successfully!")
        call message("===========================================================")


    end subroutine

    subroutine dealiasFields(this)
        class(igrid), intent(inout) :: this

        call this%spectC%dealias(this%uhat)
        call this%spectC%dealias(this%vhat)
        if (this%PeriodicInZ) then
            call transpose_y_to_z(this%what, this%cbuffzE(:,:,:,1), this%sp_gpE)
            call this%spectC%dealias_edgeField(this%cbuffzE(:,:,:,1))
            call transpose_z_to_y(this%cbuffzE(:,:,:,1),this%what,this%sp_gpE)
        else
            call this%spectE%dealias(this%what)
        end if 
        if (this%isStratified .or. this%initspinup) call this%spectC%dealias(this%That)
        if (this%UseDealiasFilterVert) then
            call this%ApplyCompactFilter()
        end if

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
        end select

    end subroutine

    subroutine reset_pointers(this, resetRHS)
        class(igrid), intent(inout), target :: this
        logical, intent(in), optional :: resetRHS 

        this%uhat => this%SfieldsC(:,:,:,1); 
        this%vhat => this%SfieldsC(:,:,:,2); 
        this%what => this%SfieldsE(:,:,:,1)
        if ((this%isStratified) .or. (this%initspinup)) then
            this%That => this%SfieldsC(:,:,:,4)
        end if

        if (present(resetRHS)) then
            if (resetRHS) then
               this%u_rhs => this%rhsC(:,:,:,1); 
               this%v_rhs => this%rhsC(:,:,:,2); 
               this%w_rhs => this%rhsE(:,:,:,1)
               if ((this%isStratified) .or. (this%initspinup)) then
                   this%T_rhs =>  this%rhsC(:,:,:,3)
               end if
            end if 
        end if 

    end subroutine

    
    subroutine TVD_RK3(this, dtforced)
        class(igrid), intent(inout), target :: this
        real(rkind), intent(in), optional :: dtforced

        if(present(dtforced)) then
          this%dt = dtforced
        else
          ! Step 0: Compute TimeStep 
          call this%compute_deltaT
        endif
       
        !!! STAGE 1
        ! First stage - everything is where it's supposed to be
        if (this%AlreadyHaveRHS) then
            this%AlreadyHaveRHS = .false.
        else
            call this%populate_rhs()
        end if
        !print*, sum(abs(this%u_rhs))
        !print*, sum(abs(this%v_rhs))
        !print*, sum(abs(this%w_rhs))
        !print*, sum(abs(this%T_rhs))
      
        this%newTimeStep = .false. 
        this%uhat1 = this%uhat + this%dt*this%u_rhs 
        this%vhat1 = this%vhat + this%dt*this%v_rhs 
        this%what1 = this%what + this%dt*this%w_rhs 
        if (this%isStratified .or. this%initspinup) this%That1 = this%That + this%dt*this%T_rhs
        ! Now set pointers so that things operate on uhat1, vhat1, etc.
        this%uhat => this%SfieldsC2(:,:,:,1); this%vhat => this%SfieldsC2(:,:,:,2); this%what => this%SfieldsE2(:,:,:,1); 
        if (this%isStratified .or. this%initspinup) this%That => this%SfieldsC2(:,:,:,3)
        ! Now perform the projection and prep for next stage
        call this%project_and_prep(this%fastCalcPressure)

        !!! STAGE 2
        ! Second stage - u, v, w are really pointing to u1, v1, w1 (which is
        ! what we want. 
        call this%populate_rhs()
        ! reset u, v, w pointers
        call this%reset_pointers()
        this%uhat1 = (3.d0/4.d0)*this%uhat + (1.d0/4.d0)*this%uhat1 + (1.d0/4.d0)*this%dt*this%u_rhs
        this%vhat1 = (3.d0/4.d0)*this%vhat + (1.d0/4.d0)*this%vhat1 + (1.d0/4.d0)*this%dt*this%v_rhs
        this%what1 = (3.d0/4.d0)*this%what + (1.d0/4.d0)*this%what1 + (1.d0/4.d0)*this%dt*this%w_rhs
        if (this%isStratified .or. this%initspinup) this%That1 = (3.d0/4.d0)*this%That + (1.d0/4.d0)*this%That1 + (1.d0/4.d0)*this%dt*this%T_rhs
        ! now set the u, v, w, pointers to u1, v1, w1
        this%uhat => this%SfieldsC2(:,:,:,1); this%vhat => this%SfieldsC2(:,:,:,2); this%what => this%SfieldsE2(:,:,:,1); 
        if (this%isStratified .or. this%initspinup) this%That => this%SfieldsC2(:,:,:,3)
        ! Now perform the projection and prep for next stage
        call this%project_and_prep(.false.)


        !!! STAGE 3 (Final Stage)
        ! Third stage - u, v, w are really pointing to u2, v2, w2 (which is what
        ! we really want. 
        call this%populate_rhs()
        ! reset u, v, w pointers
        call this%reset_pointers()
        this%uhat = (1.d0/3.d0)*this%uhat + (2.d0/3.d0)*this%uhat1 + (2.d0/3.d0)*this%dt*this%u_rhs
        this%vhat = (1.d0/3.d0)*this%vhat + (2.d0/3.d0)*this%vhat1 + (2.d0/3.d0)*this%dt*this%v_rhs
        this%what = (1.d0/3.d0)*this%what + (2.d0/3.d0)*this%what1 + (2.d0/3.d0)*this%dt*this%w_rhs
        if (this%isStratified .or. this%initspinup) this%That = (1.d0/3.d0)*this%That + (2.d0/3.d0)*this%That1 + (2.d0/3.d0)*this%dt*this%T_rhs
        ! Now perform the projection and prep for next time step
        call this%project_and_prep(.false.)


        ! Wrap up this time step 
        call this%wrapup_timestep() 
    
    end subroutine


    subroutine SSP_RK45(this, dtforced)
        class(igrid), intent(inout), target :: this
        real(rkind), intent(in), optional :: dtforced
        real(rkind), parameter :: b01 = 0.39175222657189d0 , b12 = 0.368410593050371d0, b23 = 0.25189177427169d0, b34 = 0.54497475022852d0
        real(rkind), parameter :: b35 = 0.06369246866629d0 , b45 = 0.22600748323690d0
        
        real(rkind), parameter :: a20 = 0.444370493651235d0, a21 = 0.555629506348765d0
        real(rkind), parameter :: a30 = 0.620101851488403d0, a32 = 0.379898148511597d0
        real(rkind), parameter :: a40 = 0.17807995439313d0 , a43 = 0.821920045606868d0
        real(rkind), parameter :: a52 = 0.517231671970585d0, a53 = 0.096059710526147d0, a54 = 0.386708617503269d0

        if(present(dtforced)) then
          this%dt = dtforced
        else
          ! Step 0: Compute TimeStep 
          call this%compute_deltaT
        endif
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!! STAGE 1 !!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! First stage - everything is where it's supposed to be
        if (this%AlreadyHaveRHS) then
            this%AlreadyHaveRHS = .false.
        else
            call this%populate_rhs()
        end if
        this%newTimeStep = .false.
        ! Set the pointers 
        this%uhat1 => this%uExtra(:,:,:,1); this%vhat1 => this%vExtra(:,:,:,1);
        this%what1 => this%wExtra(:,:,:,1); this%That1 => this%TExtra(:,:,:,1);
        ! Do the time step
        this%uhat1 = this%uhat + b01*this%dt*this%u_rhs 
        this%vhat1 = this%vhat + b01*this%dt*this%v_rhs 
        this%what1 = this%what + b01*this%dt*this%w_rhs 
        if (this%isStratified .or. this%initspinup) this%That1 = this%That + b01*this%dt*this%T_rhs
        ! Now set pointers so that things operate on uhat1, vhat1, etc.
        this%uhat => this%uExtra(:,:,:,1); this%vhat => this%vExtra(:,:,:,1)
        this%what => this%wExtra(:,:,:,1); this%That => this%TExtra(:,:,:,1)
        ! Now perform the projection and prep for next stage
        call this%project_and_prep(this%fastCalcPressure)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!! STAGE 2 !!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call this%populate_rhs()
        ! reset u, v, w pointers
        call this%reset_pointers()
        ! Set the new pointers
        this%uhat2 => this%uExtra(:,:,:,2); this%vhat2 => this%vExtra(:,:,:,2)
        this%what2 => this%wExtra(:,:,:,2); this%That2 => this%TExtra(:,:,:,2)
        ! Do the time step
        this%uhat2 = a20*this%uhat + a21*this%uhat1 + b12*this%dt*this%u_rhs
        this%vhat2 = a20*this%vhat + a21*this%vhat1 + b12*this%dt*this%v_rhs
        this%what2 = a20*this%what + a21*this%what1 + b12*this%dt*this%w_rhs
        if (this%isStratified .or. this%initspinup) this%That2 = a20*this%That + a21*this%That1 + b12*this%dt*this%T_rhs
        ! now set the u, v, w, pointers to u2, v2, w2
        this%uhat => this%uExtra(:,:,:,2); this%vhat => this%vExtra(:,:,:,2)
        this%what => this%wExtra(:,:,:,2); this%That => this%TExtra(:,:,:,2)
        ! Now perform the projection and prep for next stage
        call this%project_and_prep(.false.)


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!! STAGE 3 !!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call this%populate_rhs()
        ! reset u, v, w pointers
        call this%reset_pointers()
        ! Set the pointers 
        this%uhat3 => this%uExtra(:,:,:,3); this%vhat3 => this%vExtra(:,:,:,3)
        this%what3 => this%wExtra(:,:,:,3); this%That3 => this%TExtra(:,:,:,3)
        ! Do the time step 
        this%uhat3 = a30*this%uhat + a32*this%uhat2 + b23*this%dt*this%u_rhs
        this%vhat3 = a30*this%vhat + a32*this%vhat2 + b23*this%dt*this%v_rhs
        this%what3 = a30*this%what + a32*this%what2 + b23*this%dt*this%w_rhs
        if (this%isStratified .or. this%initspinup) this%That3 = a30*this%That + a32*this%That2 + b23*this%dt*this%T_rhs
        ! now set u, v, w pointers to point to u3, v3, w3
        this%uhat => this%uExtra(:,:,:,3); this%vhat => this%vExtra(:,:,:,3)
        this%what => this%wExtra(:,:,:,3); this%That => this%TExtra(:,:,:,3)
        ! Now perform the projection and prep for next time step
        call this%project_and_prep(.false.)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!! STAGE 4 !!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call this%populate_rhs()
        ! reset u, v, w pointers
        call this%reset_pointers()
        ! Set the pointers 
        this%uhat4 => this%uhat; this%vhat4 => this%vhat
        this%what4 => this%what; this%That4 => this%That
        ! Do the time step 
        this%uhat4 = a40*this%uhat + a43*this%uhat3 + b34*this%dt*this%u_rhs
        this%vhat4 = a40*this%vhat + a43*this%vhat3 + b34*this%dt*this%v_rhs
        this%what4 = a40*this%what + a43*this%what3 + b34*this%dt*this%w_rhs
        if (this%isStratified .or. this%initspinup) this%That4 = a40*this%That + a43*this%That3 + b34*this%dt*this%T_rhs
        ! now set u, v, w pointers to point to u4, v4, w4
        ! < no need to do anything here since u4 is already pointing to u > 
        ! Now perform the projection and prep for next time step
        call this%project_and_prep(.false.)


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!! STAGE 5 !!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! First reset urhs pointers to point to spare buffers
        this%u_rhs => this%uRHSExtra(:,:,:,1); this%v_rhs => this%vRHSExtra(:,:,:,1);
        this%w_rhs => this%wRHSExtra(:,:,:,1); this%T_rhs => this%TRHSExtra(:,:,:,1);
        call this%populate_rhs()
        ! Reset pointers
        call this%reset_pointers(resetRHS=.true.)
        ! Do the time step 
        this%uhat = a52*this%uhat2 + a53*this%uhat3 + b35*this%dt*this%u_rhs + a54*this%uhat + b45*this%dt*this%uRHSExtra(:,:,:,1)
        this%vhat = a52*this%vhat2 + a53*this%vhat3 + b35*this%dt*this%v_rhs + a54*this%vhat + b45*this%dt*this%vRHSExtra(:,:,:,1)
        this%what = a52*this%what2 + a53*this%what3 + b35*this%dt*this%w_rhs + a54*this%what + b45*this%dt*this%wRHSExtra(:,:,:,1)
        if (this%isStratified .or. this%initspinup) this%That = a52*this%That2 + a53*this%That3 + b35*this%dt*this%T_rhs + a54*this%That + b45*this%dt*this%TRHSExtra(:,:,:,1) 
        ! Now perform the projection and prep for next time step
        call this%project_and_prep(.false.)

        ! Wrap up this time step 
        call this%wrapup_timestep() 

    end subroutine


    subroutine computePressure(this)
        class(igrid), intent(inout) :: this
      
        ! STEP 1: Populate RHS 
        call this%populate_rhs()

        ! STEP 2: Compute pressure
        if (this%fastCalcPressure) then
            call this%Padepoiss%getPressureAndUpdateRHS(this%u_rhs,this%v_rhs,this%w_rhs,this%pressure)
        else
            call this%padepoiss%getPressure(this%u_rhs,this%v_rhs,this%w_rhs,this%pressure)
        end if 

        if (.not. useSkewSymm) then 
            ! You are using the rotational form. 
            ! This means that the pressure is really 
            ! the Bernoulli pressure. Need to subtract 
            ! out the kinetic energy. 
            call this%correctPressureRotationalForm()
        end if

        ! STEP 3: Inform the other subroutines that you already have RHS
        this%AlreadyHaveRHS = .true. 

    end subroutine

    subroutine compute_deltaT(this)
        use reductions, only: p_maxval
        class(igrid), intent(inout), target :: this
        real(rkind) :: TSmax , Tvisc
        real(rkind), dimension(:,:,:), pointer :: rb1, rb2
        
        rb1 => this%rbuffxC(:,:,:,1)
        rb2 => this%rbuffxC(:,:,:,2)

        if (this%useCFL) then
            rb1 = (one/this%dx)*this%u
            rb2 = (one/this%dy)*this%v 
            rb1 = abs(rb1) 
            rb2 = abs(rb2)
            rb1 = rb1 + rb2
            rb2 = (one/this%dz)*this%wC
            rb2 = abs(rb2)
            rb1 = rb1 + rb2
            TSmax = p_maxval(rb1)
            this%dt = this%CFL/TSmax

            if (.not. this%isInviscid) then
               Tvisc = this%CFL*this%Re*(min(this%dx,this%dy,this%dz)**2)
               this%dt = min(this%dt,Tvisc)
            end if
        end if 


    end subroutine


    function getMaxKE(this) result(maxKE)
        class(igrid), intent(inout) :: this
        real(rkind)  :: maxKE

        this%rbuffxC(:,:,:,1) = this%u**2 + this%v**2 + this%wC**2
        maxKE = half*p_maxval(maxval(this%rbuffxC))

    end function

   function getMeanKE(this) result(meanTKE)
        class(igrid), intent(inout) :: this
        real(rkind)  :: meanTKE

        this%rbuffxC(:,:,:,1) = this%u**2 + this%v**2 + this%wC**2
        meanTKE = half*p_sum(sum(this%rbuffxC))/(real(this%nx)*real(this%ny)*real(this%nz))

   end function

    subroutine interp_PrimitiveVars(this)
        class(igrid), intent(inout), target :: this
        complex(rkind), dimension(:,:,:), pointer :: ybuffC, zbuffC, zbuffE
        
        zbuffE => this%cbuffzE(:,:,:,1)
        zbuffC => this%cbuffzC(:,:,:,1)
        ybuffC => this%cbuffyC(:,:,:,1)

        ! Step 1: Interpolate w -> wC
        call transpose_y_to_z(this%what,zbuffE,this%sp_gpE)
        call this%Pade6opZ%interpz_E2C(zbuffE,zbuffC,wBC_bottom, wBC_top)
        call transpose_z_to_y(zbuffC,this%whatC,this%sp_gpC)
        call this%spectC%ifft(this%whatC,this%wC)

        ! Step 2: Interpolate u -> uE
        call transpose_y_to_z(this%uhat,zbuffC,this%sp_gpC)
        call this%Pade6opZ%interpz_C2E(zbuffC,zbuffE,uBC_bottom, uBC_top)
        call transpose_z_to_y(zbuffE,this%uEhat, this%sp_gpE)
        call this%spectE%ifft(this%uEhat, this%uE)

        ! Step 3: Interpolate v -> vE
        call transpose_y_to_z(this%vhat,zbuffC,this%sp_gpC)
        call this%Pade6opZ%interpz_C2E(zbuffC,zbuffE,vBC_bottom, vBC_top)
        call transpose_z_to_y(zbuffE,this%vEhat, this%sp_gpE)
        call this%spectE%ifft(this%vEhat, this%vE)
        

        ! Step 4: Interpolate T
        if (this%isStratified .or. this%initspinup) then
            call transpose_y_to_z(this%That,zbuffC,this%sp_gpC)
            call this%Pade6opZ%interpz_C2E(zbuffC,zbuffE,TBC_bottom, TBC_top)
            if (this%botBC_Temp == 0) then 
                zbuffE(:,:,1) = zero 
                if (nrank == 0) then
                    zbuffE(1,1,1) = this%Tsurf*real(this%nx,rkind)*real(this%ny,rkind)
                end if 
            end if
            call transpose_z_to_y(zbuffE,this%TEhat,this%sp_gpE)
            call this%spectE%ifft(this%TEhat,this%TE)
        end if 
    end subroutine

    subroutine printDivergence(this)
        class(igrid), intent(inout) :: this
        call this%padepoiss%DivergenceCheck(this%uhat, this%vhat, this%what, this%divergence)
    end subroutine 

    subroutine destroy(this)
        class(igrid), intent(inout) :: this
      
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
    end subroutine

    subroutine addNonLinearTerm_Rot(this)
        class(igrid), intent(inout), target :: this
        real(rkind),    dimension(:,:,:), pointer :: dudy, dudz, dudx
        real(rkind),    dimension(:,:,:), pointer :: dvdx, dvdy, dvdz
        real(rkind),    dimension(:,:,:), pointer :: dwdx, dwdy, dwdz
        real(rkind),    dimension(:,:,:), pointer :: dvdzC, dudzC
        real(rkind),    dimension(:,:,:), pointer :: dwdxC, dwdyC
        real(rkind),    dimension(:,:,:), pointer :: T1C, T2C, T1E, T2E 
        complex(rkind), dimension(:,:,:), pointer :: fT1C, fT2C, fT1E, fT2E 
        complex(rkind), dimension(:,:,:), pointer :: tzC, tzE

        dudx  => this%duidxjC(:,:,:,1); dudy  => this%duidxjC(:,:,:,2); dudzC => this%duidxjC(:,:,:,3); 
        dvdx  => this%duidxjC(:,:,:,4); dvdy  => this%duidxjC(:,:,:,5); dvdzC => this%duidxjC(:,:,:,6); 
        dwdxC => this%duidxjC(:,:,:,7); dwdyC => this%duidxjC(:,:,:,8); dwdz  => this%duidxjC(:,:,:,9); 

        !dwdx => this%duidxjE(:,:,:,1); dwdy => this%duidxjE(:,:,:,2);
        !dudz => this%duidxjE(:,:,:,3); dvdz => this%duidxjE(:,:,:,4);

        dwdx => this%duidxjE(:,:,:,7); dwdy => this%duidxjE(:,:,:,8);
        dudz => this%duidxjE(:,:,:,3); dvdz => this%duidxjE(:,:,:,6);

        T1C => this%rbuffxC(:,:,:,1); T2C => this%rbuffxC(:,:,:,2)
        T1E => this%rbuffxE(:,:,:,1); T2E => this%rbuffxE(:,:,:,2)
       
        fT1C => this%cbuffyC(:,:,:,1); fT2C => this%cbuffyC(:,:,:,2)
        fT1E => this%cbuffyE(:,:,:,1); fT2E => this%cbuffyE(:,:,:,2)
        
        tzC => this%cbuffzC(:,:,:,1); tzE => this%cbuffzE(:,:,:,1)


        T1C = dvdx - dudy
        T1C = T1c*this%v
        call this%spectC%fft(T1C,fT1C)
        T2E = dwdx - dudz
        T2E = T2E*this%w
        call this%spectE%fft(T2E,fT2E)
        call transpose_y_to_z(fT2E,tzE, this%sp_gpE)
        call this%Pade6opZ%interpz_E2C(tzE,tzC,0,0)
        call transpose_z_to_y(tzC,this%u_rhs, this%sp_gpC)
        this%u_rhs = this%u_rhs + fT1C


        T1C = dudy - dvdx
        T1C = T1C*this%u
        call this%spectC%fft(T1C,fT1C)
        T2E = dwdy - dvdz
        T2E = T2E*this%w
        call this%spectE%fft(T2E,fT2E)
        call transpose_y_to_z(fT2E,tzE, this%sp_gpE)
        call this%Pade6opZ%interpz_E2C(tzE,tzC,0,0)
        call transpose_z_to_y(tzC,this%v_rhs, this%sp_gpC)
        this%v_rhs = this%v_rhs + fT1C

        T1E = dudz - dwdx
        T1E = T1E*this%uE
        T2E = dvdz - dwdy
        T2E = T2E*this%vE
        T1E = T1E + T2E
        call this%spectE%fft(T1E,this%w_rhs)

        if (this%isStratified .or. this%initspinup) then
            T1C = -this%u*this%dTdxC 
            T2C = -this%v*this%dTdyC
            T1C = T1C + T2C
            call this%spectC%fft(T1c,fT1C)
            T1E = -this%w*this%dTdzE
            call this%spectE%fft(T1E,fT1E)
            call transpose_y_to_z(fT2E,tzE, this%sp_gpE)
            call this%Pade6opZ%interpz_E2C(tzE,tzC,WdTdzBC_bottom,WdTdzBC_top)
            call transpose_z_to_y(tzC,this%T_rhs, this%sp_gpC)
            this%T_rhs = this%T_rhs + fT1C
        end if

    end subroutine

    subroutine addNonLinearTerm_skewSymm(this)
        class(igrid), intent(inout), target :: this
        real(rkind),    dimension(:,:,:), pointer :: dudy, dudz, dudx
        real(rkind),    dimension(:,:,:), pointer :: dvdx, dvdy, dvdz
        real(rkind),    dimension(:,:,:), pointer :: dwdx, dwdy, dwdz
        real(rkind),    dimension(:,:,:), pointer :: dvdzC, dudzC
        real(rkind),    dimension(:,:,:), pointer :: dwdxC, dwdyC
        real(rkind),    dimension(:,:,:), pointer :: T1C, T2C, T1E, T2E 
        complex(rkind), dimension(:,:,:), pointer :: fT1C, fT2C, fT1E, fT2E 
        complex(rkind), dimension(:,:,:), pointer :: tzC, tzE

        dudx  => this%duidxjC(:,:,:,1); dudy  => this%duidxjC(:,:,:,2); dudzC => this%duidxjC(:,:,:,3); 
        dvdx  => this%duidxjC(:,:,:,4); dvdy  => this%duidxjC(:,:,:,5); dvdzC => this%duidxjC(:,:,:,6); 
        dwdxC => this%duidxjC(:,:,:,7); dwdyC => this%duidxjC(:,:,:,8); dwdz  => this%duidxjC(:,:,:,9); 

        dwdx => this%duidxjE(:,:,:,7); dwdy => this%duidxjE(:,:,:,8);
        dudz => this%duidxjE(:,:,:,3); dvdz => this%duidxjE(:,:,:,6);

        T1C => this%rbuffxC(:,:,:,1); T2C => this%rbuffxC(:,:,:,2)
        T1E => this%rbuffxE(:,:,:,1); T2E => this%rbuffxE(:,:,:,2)
       
        fT1C => this%cbuffyC(:,:,:,1); fT2C => this%cbuffyC(:,:,:,2)
        fT1E => this%cbuffyE(:,:,:,1); fT2E => this%cbuffyE(:,:,:,2)
        
        tzC => this%cbuffzC(:,:,:,1); tzE => this%cbuffzE(:,:,:,1)


        T1C = dudx*this%u
        T2C = dudy*this%v
        T1C = T1C + T2C
        T1E = dudz*this%w
        call this%spectC%fft(T1C,fT1C)
        call this%spectE%fft(T1E,fT1E)
        call transpose_y_to_z(fT1E,tzE, this%sp_gpE)
        call this%Pade6opZ%interpz_E2C(tzE,tzC,WdUdzBC_bottom,WdUdzBC_top)
        call transpose_z_to_y(tzC,this%u_rhs, this%sp_gpC)
        this%u_rhs = this%u_rhs + fT1C
        
        T1C = dvdx*this%u
        T2C = dvdy*this%v
        T1C = T1C + T2C
        T1E = dvdz*this%w
        call this%spectC%fft(T1C,fT1C)
        call this%spectE%fft(T1E,fT1E)
        call transpose_y_to_z(fT1E,tzE, this%sp_gpE)
        call this%Pade6opZ%interpz_E2C(tzE,tzC,WdVdzBC_bottom,WdVdzBC_top)
        call transpose_z_to_y(tzC,this%v_rhs, this%sp_gpC)
        this%v_rhs = this%v_rhs + fT1C
        
        T1E = dwdx*this%uE
        T2E = dwdy*this%vE
        T2E = T1E + T2E
        call this%spectE%fft(T2E,fT2E)
        T1C = dwdz*this%wC
        call this%spectC%fft(T1C,fT1C)
        call transpose_y_to_z(fT1C,tzC, this%sp_gpC)
        call this%Pade6opZ%interpz_C2E(tzC,tzE,WdWdzBC_bottom,WdWdzBC_top)
        call transpose_z_to_y(tzE,this%w_rhs, this%sp_gpE)
        this%w_rhs = this%w_rhs + fT2E

        T1C = this%u*this%u
        call this%spectC%fft(T1C,fT1C)
        call this%spectC%mtimes_ik1_ip(fT1C)
        this%u_rhs = this%u_rhs + fT1C

        T1C = this%v*this%v
        call this%spectC%fft(T1C,fT1C)
        call this%spectC%mtimes_ik2_ip(fT1C)
        this%v_rhs = this%v_rhs + fT1C

        T1C = this%wC*this%wC
        call this%spectC%fft(T1C,fT1C)
        call transpose_y_to_z(fT1C,tzC,this%sp_gpC)
        call this%Pade6opZ%ddz_C2E(tzC,tzE,WWBC_bottom,WWBC_top)
        call transpose_z_to_y(tzE,fT1E,this%sp_gpE)
        this%w_rhs = this%w_rhs + fT1E

        T1C = this%u*this%v
        call this%spectC%fft(T1C,fT1C)
        call this%spectC%mtimes_ik2_oop(fT1C,fT2C)
        this%u_rhs = this%u_rhs + fT2C
        call this%spectC%mtimes_ik1_ip(fT1C)
        this%v_rhs = this%v_rhs + fT1C

        T1E = this%uE*this%w
        call this%spectE%fft(T1E,fT1E)
        call transpose_y_to_z(fT1E,TzE,this%sp_gpE)
        call this%Pade6opZ%ddz_E2C(tzE,tzC,UWBC_bottom,UWBC_top)
        call transpose_z_to_y(tzC,fT1C,this%sp_gpC)
        this%u_rhs = this%u_rhs + fT1C
        
        call this%spectE%mtimes_ik1_ip(fT1E)
        this%w_rhs = this%w_rhs + fT1E


        T1E = this%vE*this%w
        call this%spectE%fft(T1E,fT1E)
        call transpose_y_to_z(fT1E,TzE,this%sp_gpE)
        call this%Pade6opZ%ddz_E2C(tzE,tzC,VWBC_bottom,VWBC_top)
        call transpose_z_to_y(tzC,fT1C,this%sp_gpC)
        this%v_rhs = this%v_rhs + fT1C

        call this%spectE%mtimes_ik2_ip(fT1E)
        this%w_rhs = this%w_rhs + fT1E

        this%u_rhs = -half*this%u_rhs
        this%v_rhs = -half*this%v_rhs
        this%w_rhs = -half*this%w_rhs


        if (this%isStratified .or. this%initspinup) then
            T1C = -this%u*this%T
            call this%spectC%fft(T1C,this%T_rhs)
            call this%spectC%mtimes_ik1_ip(this%T_rhs)
            T1C = -this%v*this%T
            call this%spectC%fft(T1C,fT1C)
            call this%spectC%mtimes_ik2_ip(fT1C)
            this%T_rhs = this%T_rhs + fT1C
            T1E = -this%w * this%TE
            call this%spectE%fft(T1E,fT1E)
            call transpose_y_to_z(fT1E,TzE,this%sp_gpE)
            call this%Pade6opZ%ddz_E2C(tzE,tzC,WTBC_bottom,WTBC_top)
            call transpose_z_to_y(tzC,fT1C,this%sp_gpC)
            this%T_rhs = this%T_rhs + fT1C
        end if 

    end subroutine

    subroutine addCoriolisTerm(this)
        class(igrid), intent(inout), target :: this
        complex(rkind), dimension(:,:,:), pointer :: ybuffE, ybuffC1, ybuffC2, zbuffC, zbuffE

        ybuffE => this%cbuffyE(:,:,:,1)
        ybuffC1 => this%cbuffyC(:,:,:,1)
        ybuffC2 => this%cbuffyC(:,:,:,2)
        zbuffE => this%cbuffzE(:,:,:,1)
        zbuffC => this%cbuffzC(:,:,:,1)


        ! u equation 
        ybuffC1    = (two/this%Ro)*(-this%coriolis_omegaY*this%whatC - this%coriolis_omegaZ*(this%GyHat - this%vhat))
        this%u_rhs = this%u_rhs  + ybuffC1
        
        ! v equation
        !ybuffC2    = (two/this%Ro)*(this%coriolis_omegaZ*(this%GxHat - this%uhat))
        ybuffC2    = (two/this%Ro)*(this%coriolis_omegaZ*(this%GxHat - this%uhat) + this%coriolis_omegaX*this%whatC)
        this%v_rhs = this%v_rhs + ybuffC2 
        
        ! w equation 
        ! The real equation is given as:
        ! this%w_rhs = this%w_rhs - this%coriolis_omegaY*(this%GxHat - this%uhat)/this%Ro
        ! But we evaluate this term as:
      
        ybuffE = (two/this%Ro)*(-this%coriolis_omegaY*(this%Gxhat_Edge - this%uEhat) + this%coriolis_omegaY*(this%Gyhat_Edge - this%vEhat)) 
        if (this%spectE%CarryingZeroK) then
            ybuffE(1,1,:) = cmplx(zero,zero,rkind)
        end if 
        this%w_rhs = this%w_rhs + ybuffE 
        ! The residual quantity (Gx - <u>)*cos(alpha)/Ro is accomodated in
        ! pressure

        if (this%storeFbody) then
            call this%spectC%ifft(ybuffC1,this%fbody_x)
            call this%spectC%ifft(ybuffC2,this%fbody_y)
            call this%spectE%ifft(ybuffE ,this%fbody_z)
        end if
    end subroutine  

    subroutine addSponge(this)
        class(igrid), intent(inout), target :: this
        complex(rkind), dimension(:,:,:), pointer :: deviationC

        deviationC => this%cbuffyC(:,:,:,1)
        
        deviationC = this%uhat - this%ubase
        this%u_rhs = this%u_rhs - (this%RdampC/this%dt)*deviationC

        deviationC = this%vhat - this%vbase
        this%v_rhs = this%v_rhs - (this%RdampC/this%dt)*deviationC
        
        this%w_rhs = this%w_rhs - (this%RdampE/this%dt)*this%what ! base value for w is zero  

        !deviationC = this%That - this%Tbase
        !this%T_rhs = this%T_rhs - (this%RdampC/this%dt)*deviationC

    end subroutine

    subroutine addExtraForcingTerm(this)
        class(igrid), intent(inout) :: this
        this%u_rhs = this%u_rhs + this%dpF_dxhat
        if (this%storeFbody) then
            this%fbody_x = this%fbody_x + this%dpFdx 
        end if
    end subroutine

    subroutine AddBuoyancyTerm(this)
        class(igrid), intent(inout), target :: this
        complex(rkind), dimension(:,:,:), pointer :: fT1E 
        real(rkind), dimension(:,:,:), pointer :: rbuffE

        fT1E => this%cbuffyE(:,:,:,1)
        rbuffE => this%rbuffxE(:,:,:,1)

        !fT1E = (this%TEhat)/(this%ThetaRef*this%Fr*this%Fr)
        fT1E = (this%TEhat)*this%BuoyancyFact ! See definition of buoyancy factor in init 
        if (this%spectE%carryingZeroK) then
            fT1E(1,1,:) = cmplx(zero,zero,rkind)
        end if 
        this%w_rhs = this%w_rhs + fT1E 
        
        if (this%useSponge) then
            call this%addSponge
        end if

        if (this%storeFbody) then
            call this%spectE%ifft(fT1E, rbuffE)
            this%fbody_z = this%fbody_z + rbuffE
        end if

    end subroutine

    subroutine populate_rhs(this)
        class(igrid), intent(inout) :: this
        !integer,           intent(in)    :: RKstage

        ! Step 1: Non Linear Term 
        if (useSkewSymm) then
            call this%addNonLinearTerm_skewSymm()
        else
            call this%AddNonLinearTerm_Rot()
        end if

        !print*, "1:" 
        !print*, "urhs:", sum(abs(this%u_rhs))
        !print*, "vrhs:", sum(abs(this%v_rhs))
        !print*, "wrhs:", sum(abs(this%w_rhs))
        !print*, "Trhs:", sum(abs(this%T_rhs))

        ! Step 2: Coriolis Term
        if (this%useCoriolis) then
            call this%AddCoriolisTerm()
        end if 
        
        !print*, "2:" 
        !print*, "urhs:", sum(abs(this%u_rhs))
        !print*, "vrhs:", sum(abs(this%v_rhs))
        !print*, "wrhs:", sum(abs(this%w_rhs))
        !print*, "Trhs:", sum(abs(this%T_rhs))
        
        ! Step 3a: Extra Forcing 
        if (this%useExtraForcing) then
            call this%addExtraForcingTerm()
        end if

        ! Step 3b: Wind Turbines
        !if (this%useWindTurbines .and. (RKstage==1)) then
        if (this%useWindTurbines) then
           if (allocated(this%inst_horz_avg_turb)) then
               call this%WindTurbineArr%getForceRHS(this%dt, this%u, this%v, this%wC,&
                                    this%u_rhs, this%v_rhs, this%w_rhs, this%newTimestep, this%inst_horz_avg_turb)
           else
               call this%WindTurbineArr%getForceRHS(this%dt, this%u, this%v, this%wC,&
                                    this%u_rhs, this%v_rhs, this%w_rhs, this%newTimestep)
           end if
        end if 
       
        !print*, "3:" 
        !print*, "urhs:", sum(abs(this%u_rhs))
        !print*, "vrhs:", sum(abs(this%v_rhs))
        !print*, "wrhs:", sum(abs(this%w_rhs))
        !!print*, "Trhs:", sum(abs(this%T_rhs))
        !print '(A,ES26.16)', "Trhs:", sum(abs(this%T_rhs))

        ! Step 4: Buoyance + Sponge (inside Buoyancy)
        if (this%isStratified .or. this%initspinup) then
            call this%addBuoyancyTerm()
        end if 
        
        !print*, "4:" 
        !print*, "urhs:", sum(abs(this%u_rhs))
        !print*, "vrhs:", sum(abs(this%v_rhs))
        !print*, "wrhs:", sum(abs(this%w_rhs))
        !!print*, "Trhs:", sum(abs(this%T_rhs))
        !print '(A,ES26.16)', "Trhs:", sum(abs(this%T_rhs))      

        ! Step 5: Viscous Term (only if simulation if NOT inviscid)
        if (.not. this%isInviscid) then
            call this%addViscousTerm()
        end if
        
        !print*, "5:" 
        !print*, "urhs:", sum(abs(this%u_rhs))
        !print*, "vrhs:", sum(abs(this%v_rhs))
        !print*, "wrhs:", sum(abs(this%w_rhs))
        !print*, "Trhs:", sum(abs(this%T_rhs))

        ! Step 6: SGS Viscous Term
        if (this%useSGS) then
            call this%sgsmodel%getRHS_SGS(this%u_rhs,      this%v_rhs, this%w_rhs,      this%duidxjC, this%duidxjE, &
                                          this%duidxjEhat, this%uEhat, this%vEhat,      this%what,    this%uhat,    &
                                          this%vhat,       this%That,  this%u,          this%v,       this%uE,      &
                                          this%vE,         this%w,     this%newTimeStep                             )

            if (this%isStratified .or. this%initspinup) then
               call this%sgsmodel%getRHS_SGS_Scalar(this%T_rhs, this%dTdxC, this%dTdyC, this%dTdzE)
            end if
            
        end if
        
        !print*, "6:" 
        !print*, "urhs:", sum(abs(this%u_rhs))
        !print*, "vrhs:", sum(abs(this%v_rhs))
        !print*, "wrhs:", sum(abs(this%w_rhs))
        !print*, "Trhs:", sum(abs(this%T_rhs))

        ! Step 7: Fringe source term if fringe is being used (non-periodic)
        if (this%usedoublefringex) then
                call this%fringe_x1%addFringeRHS(this%dt, this%u_rhs, this%v_rhs, this%w_rhs, this%u, this%v, this%w)
                call this%fringe_x2%addFringeRHS(this%dt, this%u_rhs, this%v_rhs, this%w_rhs, this%u, this%v, this%w)
        else
            if (this%useFringe) then
                call this%fringe_x%addFringeRHS(this%dt, this%u_rhs, this%v_rhs, this%w_rhs, this%u, this%v, this%w)
            end if 
        end if 

        ! Step 8: HIT forcing source term
        if (this%useHITForcing) then
            call this%hitforce%getRHS_HITForcing(this%u_rhs, this%v_rhs, this%w_rhs, this%uhat, this%vhat, this%what, this%newTimeStep)
        end if 

        !if (nrank == 0) print*, maxval(abs(this%u_rhs)), maxval(abs(this%v_rhs)), maxval(abs(this%w_rhs))
    end subroutine

    subroutine addViscousTerm(this)
        class(igrid), intent(inout) :: this
        integer :: i, j, k
        real(rkind) :: oneByRe, molecularDiff
        complex(rkind) :: tmp1, tmp2

        oneByRe = one/this%Re

        do k = 1,size(this%u_rhs,3)
           do j = 1,size(this%u_rhs,2)
              !$omp simd
              do i = 1,size(this%u_rhs,1)
                  tmp1 = -this%spectC%kabs_sq(i,j,k)*this%uhat(i,j,k) + this%d2udz2hatC(i,j,k)
                  tmp2 = -this%spectC%kabs_sq(i,j,k)*this%vhat(i,j,k) + this%d2vdz2hatC(i,j,k)
                  this%u_rhs(i,j,k) = this%u_rhs(i,j,k) + oneByRe*tmp1
                  this%v_rhs(i,j,k) = this%v_rhs(i,j,k) + oneByRe*tmp2
               end do
            end do
         end do

        do k = 1,size(this%w_rhs,3)
           do j = 1,size(this%w_rhs,2)
              !$omp simd
              do i = 1,size(this%w_rhs,1)
                  tmp1 = -this%spectE%kabs_sq(i,j,k)*this%what(i,j,k) + this%d2wdz2hatE(i,j,k)
                  this%w_rhs(i,j,k) = this%w_rhs(i,j,k) + oneByRe*tmp1
               end do
            end do
         end do

         if (this%isStratified) then
            molecularDiff = one/(this%Re*this%PrandtlFluid)
            do k = 1,size(this%T_rhs,3)
               do j = 1,size(this%T_rhs,2)
                  !$omp simd
                  do i = 1,size(this%T_rhs,1)
                     tmp1 = -this%spectC%kabs_sq(i,j,k)*this%That(i,j,k) + this%d2Tdz2hatC(i,j,k) 
                     this%T_rhs(i,j,k) = this%T_rhs(i,j,k) + molecularDiff*tmp1
                  end do 
               end do 
            end do 
            
         end if
        !this%cbuffyC(:,:,:,1) = -this%spectC%kabs_sq*this%uhat + this%d2udz2hatC
        !this%u_rhs = this%u_rhs + (one/this%Re)*this%cbuffyC(:,:,:,1)

        !this%cbuffyC(:,:,:,1) = -this%spectC%kabs_sq*this%vhat + this%d2vdz2hatC
        !this%v_rhs = this%v_rhs + (one/this%Re)*this%cbuffyC(:,:,:,1)

        !this%cbuffyE(:,:,:,1) = -this%spectE%kabs_sq*this%what + this%d2wdz2hatE
        !this%w_rhs = this%w_rhs + (one/this%Re)*this%cbuffyE(:,:,:,1)

    end subroutine


    subroutine project_and_prep(this, AlreadyProjected)
        class(igrid), intent(inout) :: this
        logical, intent(in) :: AlreadyProjected

        ! Step 1: Dealias
        !call this%spectC%dealias(this%uhat)
        !call this%spectC%dealias(this%vhat)
        !call this%spectE%dealias(this%what)
        !if (this%isStratified .or. this%initspinup) call this%spectC%dealias(this%That)
        !if (this%UseDealiasFilterVert) then
        !    call this%ApplyCompactFilter()
        !end if
        call this%dealiasFields()

        ! Step 2: Pressure projection
        if (.not. AlreadyProjected) then
            call this%padepoiss%PressureProjection(this%uhat,this%vhat,this%what)
            if (mod(this%step,this%t_DivergenceCheck) == 0) then
                call this%padepoiss%DivergenceCheck(this%uhat, this%vhat, this%what, this%divergence,.true.)
            end if 
        end if 

        ! Step 3: Take it back to physical fields
        call this%spectC%ifft(this%uhat,this%u)
        call this%spectC%ifft(this%vhat,this%v)
        call this%spectE%ifft(this%what,this%w)
        if (this%isStratified .or. this%initspinup) call this%spectC%ifft(this%That,this%T)
    
        ! STEP 4: Interpolate the cell center values of w
        !call this%compute_and_bcast_surface_Mn()
        call this%interp_PrimitiveVars()

        ! STEP 5: Compute duidxjC 
        call this%compute_duidxj()
        if (this%isStratified .or. this%initspinup) call this%compute_dTdxi() 

    end subroutine


    subroutine updateProbes(this)
        class(igrid), intent(inout) :: this
        integer :: idx

        if (this%doIhaveAnyProbes) then
            do idx = 1,this%nprobes
                this%probe_data(1,idx,this%step) = this%tsim
                this%probe_data(2,idx,this%step) = this%u (this%probes(1,idx),this%probes(2,idx),this%probes(3,idx))
                this%probe_data(3,idx,this%step) = this%v (this%probes(1,idx),this%probes(2,idx),this%probes(3,idx))
                this%probe_data(4,idx,this%step) = this%wC(this%probes(1,idx),this%probes(2,idx),this%probes(3,idx))
                if (this%isStratified) then
                    this%probe_data(5,idx,this%step) = this%T(this%probes(1,idx),this%probes(2,idx),this%probes(3,idx))
                end if
                if (this%fastCalcPressure) then
                    this%probe_data(6,idx,this%step) = this%Pressure(this%probes(1,idx),this%probes(2,idx),this%probes(3,idx))
                end if
            end do 
        end if

        ! KS - preprocess
        if (this%PreprocessForKS) then
            if (.not. this%KSupdated) then
                call this%LES2KS%applyFilterForKS(this%u, this%v, this%wC)
                this%KSupdated = .true. 
            end if
            if (this%doIhaveAnyProbes) then
                do idx = 1,this%nprobes
                    this%KS_probe_data(1,idx,this%step) = this%tsim
                    this%KS_probe_data(2,idx,this%step) = this%ufil4KS (this%probes(1,idx),this%probes(2,idx),this%probes(3,idx))
                    this%KS_probe_data(3,idx,this%step) = this%vfil4KS (this%probes(1,idx),this%probes(2,idx),this%probes(3,idx))
                    this%KS_probe_data(4,idx,this%step) = this%wfil4KS(this%probes(1,idx),this%probes(2,idx),this%probes(3,idx))
                end do 
            end if
        end if

    end subroutine


    subroutine dumpProbes(this)
        use basic_io, only: write_2d_ascii
        class(igrid), intent(in) :: this
        character(len=clen) :: tempname, fname
        integer :: pid, idx

        do idx = 1,this%nprobes
            pid = this%probes(4,idx)
            write(tempname,"(A3,I2.2,A6,I3.3,A4,I6.6,A4,I6.6,A4)") "Run",this%runID, "_PROBE",pid,"_tst",this%probeStartStep,"_ten",this%step,".out"
            fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
            call write_2d_ascii(transpose(this%probe_data(:,idx,this%probeStartStep:this%step)), fname)
        end do 

        ! KS - preprocess
        if (this%PreprocessForKS) then
            do idx = 1,this%nprobes
                pid = this%probes(4,idx)
                write(tempname,"(A3,I2.2,A9,I3.3,A4,I6.6,A4,I6.6,A4)") "Run",this%runID, "_PROBE_KS",pid,"_tst",this%probeStartStep,"_ten",this%step,".out"
                fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                call write_2d_ascii(transpose(this%KS_probe_data(:,idx,this%probeStartStep:this%step)), fname)
            end do 
        end if
        
    end subroutine

    subroutine wrapup_timestep(this)
        class(igrid), intent(inout) :: this

        logical :: forceWrite, exitStat, forceDumpPressure, restartWrite, forceDumpProbes 
        integer :: ierr = -1, ierr2

        ! STEP 1: Update Time, BCs and record probe data
        this%step = this%step + 1; this%tsim = this%tsim + this%dt
        this%newTimeStep = .true. 
        if (this%isStratified) then
            if (this%botBC_Temp == 0) then
                this%Tsurf = this%Tsurf0 + this%dTsurf_dt*this%tsim
            end if
        end if  
        if (this%PreprocessForKS) this%KSupdated = .false. 
        if (this%useProbes) call this%updateProbes()
        if (this%computevorticity) call this%compute_vorticity

        ierr = -1; forceWrite = .FALSE.; exitstat = .FALSE.; forceDumpPressure = .FALSE.; 
        forceDumpProbes = .false.; restartWrite = .FALSE. 

        if(this%tsim > this%tstop) then
          forceWrite = .TRUE.
          restartWrite = .TRUE.
          if (this%useProbes) forceDumpProbes = .TRUE.
          call message(0,"The simulation has ended.")
          call message(1,"Dumping a restart file.")
        endif

        if (this%useSystemInteractions) then
            if (mod(this%step,this%tSystemInteractions) == 0) then
                exitStat = check_exit(this%controlDir)
                if (exitStat) forceWrite = .true.
                open(777,file=trim(this%controlDir)//"/dumppdo",status='old',iostat=ierr)
                if(ierr==0) then
                    forceWrite = .TRUE.
                    call message(1, "Forced Dump because found file dumppdo")
                    call message(2, "Current Time Step is:", this%step)
                    if (nrank .ne. 0) close(777)
                    call mpi_barrier(mpi_comm_world, ierr2)
                    if(nrank==0) close(777, status='delete')
                else
                    close(777)
                endif
           

                if (this%useProbes) then
                    open(777,file=trim(this%controlDir)//"/dumpprobes",status='old',iostat=ierr)
                    if (ierr==0) then
                        forceDumpProbes = .true. 
                        call message(1, "Forced Dump for PROBES because found file dumpprobes")
                        call message(2, "Current Time Step is:", this%step)
                        if (nrank .ne. 0) close(777)
                        call mpi_barrier(mpi_comm_world, ierr2)
                        if(nrank==0) close(777, status='delete')
                    else
                        close(777)
                    end if
                end if

                if (this%storePressure) then 
                    open(777,file=trim(this%controlDir)//"/prsspdo",status='old',iostat=ierr)
                    if (ierr == 0) then
                        forceDumpPressure = .true.
                        call message(1,"Force Dump for pressure reqested; file prsspdo found.")
                        if (nrank .ne. 0) close(777)
                        call mpi_barrier(mpi_comm_world, ierr2)
                        if (nrank==0) close(777, status='delete')
                    else
                        close(777)
                    end if
                end if

                open(777,file=trim(this%controlDir)//"/dumprestart",status='old',iostat=ierr)
                if(ierr==0) then
                    restartWrite = .TRUE.
                    call message(1, "Restart Dump because found file dumprestart")
                    call message(2, "Current Time Step is:", this%step)
                    if (nrank .ne. 0) close(777)
                    call mpi_barrier(mpi_comm_world, ierr2)
                    if(nrank==0) close(777, status='delete')
                else
                    close(777)
                endif

            end if 
        end if 

        ! STEP 2: Do logistical stuff
        if (this%fastCalcPressure) then
            call this%computePressure()
        else
            if ((this%storePressure)) then
                if ((mod(this%step,this%P_compFreq)==0) .or. (forceDumpPressure)) then
                    call this%computePressure()
                end if 
                if ( (mod(this%step,this%P_dumpFreq) == 0).or. (forceDumpPressure)) then
                    call this%dumpFullField(this%pressure,"prss")
                end if 
            end if 
        end if

       if ((forceWrite.or.(mod(this%step,this%tid_compStats)==0)).and.(this%tsim > this%tSimStartStats) ) then
           if (this%timeAvgFullFields) then
               ! rhs needs to be evaluated before computing statistics to ensure
               ! that tauSGS are consistent with the velocities and velocity
               ! derivatives, which is needed for correct SGS dissipation
               if (.not. this%AlreadyHaveRHS) then
                   call this%populate_rhs()
                   this%AlreadyHaveRHS = .true. 
               end if 

               call this%compute_stats3D()
           else
           !    call this%compute_stats()
           end if 
       end if 

       if ((forceWrite.or.(mod(this%step,this%tid_statsDump)==0)).and.(this%tsim > this%tSimStartStats) ) then
           if (this%timeAvgFullFields) then
                call this%dump_stats3D()
                call mpi_barrier(mpi_comm_world, ierr)
                !stop
            else
                !call this%compute_stats()
                !call this%dump_stats()
            end if
        end if 

        if ( restartWrite .or. (mod(this%step,this%t_restartDump) == 0) ) then
            call this%dumpRestartfile()
            call message(0,"Scheduled restart file dumped.")
        end if
        
        if ( (forceWrite .or. ((mod(this%step,this%t_planeDump) == 0) .and. &
                 (this%step .ge. this%t_start_planeDump) .and. (this%step .le. this%t_stop_planeDump))) .and. (this%dumpPlanes)) then
            if (this%PreprocessForKS) then
                if (.not. this%KSupdated) then
                    call this%LES2KS%applyFilterForKS(this%u, this%v, this%wC)
                    this%KSupdated = .true. 
                end if
            end if
            call this%dump_planes()
        end if 

        !if ( (forceWrite .or. (mod(this%step,this%t_dumpKSprep) == 0)) .and. this%PreprocessForKS ) then
        !    call this%LES2KS%LES_TO_KS(this%uE,this%vE,this%w,this%step)
        !    call this%LES2KS%LES_FOR_KS(this%uE,this%vE,this%w,this%step)
        !end if 

        ! ADITYA -> NIRANJAN: You need to fix the dump_pointProbes call. For
        ! some reason, it seems to segfault for some problems. 
        !if ( (forceWrite .or. ((mod(this%step,this%t_pointProbe) == 0) .and. &
        !         (this%step .ge. this%t_start_pointProbe) .and. (this%step .le. this%t_stop_pointProbe))) .and. (this%t_pointProbe > 0)) then
        !    call this%dump_pointProbes()
        !end if 

        if (mod(this%step,this%t_dataDump) == 0) then
           call message(0,"Scheduled visualization dump.")
           call this%dumpFullField(this%u,'uVel')
           call this%dumpFullField(this%v,'vVel')
           call this%dumpFullField(this%wC,'wVel')
           call this%dumpVisualizationInfo()
           if (this%isStratified .or. this%initspinup) call this%dumpFullField(this%T,'potT')
           if (this%fastCalcPressure) call this%dumpFullField(this%pressure,'prss')
           if (this%computevorticity) then
               call this%dumpFullField(this%ox,'omgX')
               call this%dumpFullField(this%oy,'omgY')
               call this%dumpFullField(this%oz,'omgZ')
           end if 
           
           ! ADITYA -> NIRANJAN : Did you mean to use this for debugging, or do
           ! you need it to be here? If you do, then it needs to go inside
           ! turbArr.  
           !if (this%useWindTurbines) then
           !    this%WindTurbineArr%dumpTurbField = .true. ! forces will be printed out at the next time step
           !    this%WindTurbineArr%step = this%step-1
           !endif

           if (this%PreProcessForKS) then
                if (.not. this%KSupdated) then
                    call this%LES2KS%applyFilterForKS(this%u, this%v, this%wC)
                    this%KSupdated = .true. 
                end if
                call this%dumpFullField(this%uFil4KS,'uFks')
                call this%dumpFullField(this%vFil4KS,'vFks')
                call this%dumpFullField(this%wFil4KS,'wFks')
           end if
           if (this%useProbes) then
                call this%dumpProbes()    
                call message(0,"Performed a scheduled dump for probes.")
           end if
        end if

        if (this%useProbes) then
            if (forceDumpProbes) then
                call this%dumpProbes()    
                call message(0,"Performed a forced dump for probes.")
            end if
        end if

        if (forceWrite) then
           call message(2,"Performing a forced visualization dump.")
           call this%dumpFullField(this%u,'uVel')
           call this%dumpFullField(this%v,'vVel')
           call this%dumpFullField(this%wC,'wVel')
           call this%dumpVisualizationInfo()
           if (this%isStratified .or. this%initspinup) call this%dumpFullField(this%T,'potT')
           if (this%fastCalcPressure) call this%dumpFullField(this%pressure,'prss')
           if (this%computevorticity) then
               call this%dumpFullField(this%ox,'omgX')
               call this%dumpFullField(this%oy,'omgY')
               call this%dumpFullField(this%oz,'omgZ')
           end if 
           
           ! ADITYA -> NIRANJAN : Did you mean to use this for debugging, or do
           ! you need it to be here? 
           !if (this%useWindTurbines) then
           !    this%WindTurbineArr%dumpTurbField = .true. ! forces will be printed out at the next time step
           !    this%WindTurbineArr%step = this%step-1
           !endif
           if (this%PreProcessForKS) then
                if (.not. this%KSupdated) then
                    call this%LES2KS%applyFilterForKS(this%u, this%v, this%wC)
                    this%KSupdated = .true. 
                end if
                call this%dumpFullField(this%uFil4KS,'uFks')
                call this%dumpFullField(this%vFil4KS,'vFks')
                call this%dumpFullField(this%wFil4KS,'wFks')
           end if
           !call output_tecplot(gp)
        end if

        if (this%initspinup) then
         if (this%tsim > this%Tstop_initspinup) then
             this%initspinup = .false.
             call message(0,"Initialization spin up turned off. No active scalar in the problem.")
         end if 
        end if 
        ! Exit if the exitpdo file was detected earlier at the beginning of this
        ! subroutine
        if(exitstat) call GracefulExit("Found exitpdo file in control directory",1234)

    end subroutine

    subroutine AdamsBashforth(this)
        class(igrid), intent(inout) :: this
        real(rkind) :: abf1, abf2

        ! Step 0: Compute TimeStep 
        call this%compute_deltaT
        this%dtRat = this%dt/this%dtOld

        ! Step 1: Get the RHS
        if (this%AlreadyHaveRHS) then
            this%AlreadyHaveRHS = .false.
        else
            call this%populate_rhs()
        end if 

        ! Step 2: Time Advance
        if (this%step == 0) then
            this%uhat = this%uhat + this%dt*this%u_rhs 
            this%vhat = this%vhat + this%dt*this%v_rhs 
            this%what = this%what + this%dt*this%w_rhs 
            if (this%isStratified .or. this%initspinup) then
                this%That = this%That + this%dt*this%T_rhs
            end if
        else
            abf1 = (one + half*this%dtRat)*this%dt
            abf2 = -half*this%dtRat*this%dt
            this%uhat = this%uhat + abf1*this%u_rhs + abf2*this%u_Orhs
            this%vhat = this%vhat + abf1*this%v_rhs + abf2*this%v_Orhs
            this%what = this%what + abf1*this%w_rhs + abf2*this%w_Orhs
            if (this%isStratified .or. this%initspinup) then
                this%That = this%That + abf1*this%T_rhs + abf2*this%T_Orhs
            end if 
        end if 

        ! Step 3: Pressure Project and prep for the next step
        call this%project_and_prep(.false.)

        ! Step 4: Store the RHS values for the next use
        this%u_Orhs = this%u_rhs; this%v_Orhs = this%v_rhs; this%w_Orhs = this%w_rhs
        if (this%isStratified .or. this%initspinup) this%T_Orhs = this%T_rhs
        this%dtOld = this%dt

        ! Step 5: Do end of time step operations (I/O, stats, etc.)
        call this%wrapup_timestep()
    end subroutine
    
    subroutine debug(this)
        class(igrid), intent(inout), target :: this
        real(rkind),    dimension(:,:,:), pointer :: dudx, dudy, dudz
        real(rkind),    dimension(:,:,:), pointer :: dvdx, dvdy, dvdz
        real(rkind),    dimension(:,:,:), pointer :: dwdx, dwdy, dwdz
        real(rkind),    dimension(:,:,:), pointer :: dvdzC, dudzC
        real(rkind),    dimension(:,:,:), pointer :: dwdxC, dwdyC
        real(rkind), dimension(:,:,:), pointer :: rbuff
        complex(rkind), dimension(:,:,:), pointer :: cbuff, dvdzH

        dudx  => this%duidxjC(:,:,:,1); dudy  => this%duidxjC(:,:,:,2); dudzC => this%duidxjC(:,:,:,3); 
        dvdx  => this%duidxjC(:,:,:,4); dvdy  => this%duidxjC(:,:,:,5); dvdzC => this%duidxjC(:,:,:,6); 
        dwdxC => this%duidxjC(:,:,:,7); dwdyC => this%duidxjC(:,:,:,8); dwdz  => this%duidxjC(:,:,:,9); 

        !dwdx => this%duidxjE(:,:,:,1); dwdy => this%duidxjE(:,:,:,2);
        !dudz => this%duidxjE(:,:,:,3); dvdz => this%duidxjE(:,:,:,4);
        
        dwdx => this%duidxjE(:,:,:,7); dwdy => this%duidxjE(:,:,:,8);
        dudz => this%duidxjE(:,:,:,3); dvdz => this%duidxjE(:,:,:,6);

        rbuff => this%rbuffxC(:,:,:,1); cbuff => this%cbuffyC(:,:,:,1)
        dvdzH => this%duidxjChat(:,:,:,6) 

        !print*, this%dt
        !call this%spectC%ifft(this%uhat,rbuff)
        !print*, rbuff(5,5,4)
        !call this%spectC%ifft(this%vhat,rbuff)
        !print*, rbuff(5,5,4)
        
        call this%spectC%ifft(this%u_rhs,rbuff)
        print*, rbuff(5,5,4)
        call this%spectC%ifft(this%v_rhs,rbuff)
        print*, rbuff(5,5,4)
    end subroutine 

    subroutine ApplyCompactFilter(this)
        class(igrid), intent(inout), target :: this
        complex(rkind), dimension(:,:,:), pointer :: zbuff1, zbuff2, zbuff3, zbuff4
        zbuff1 => this%cbuffzC(:,:,:,1)
        zbuff2 => this%cbuffzC(:,:,:,2)
        zbuff3 => this%cbuffzE(:,:,:,1)
        zbuff4 => this%cbuffzE(:,:,:,2)

        call transpose_y_to_z(this%uhat,zbuff1, this%sp_gpC)
        call this%filzC%filter3(zbuff1,zbuff2,this%nxZ, this%nyZ)
        call transpose_z_to_y(zbuff1,this%uhat, this%sp_gpC)

        call transpose_y_to_z(this%vhat,zbuff1, this%sp_gpC)
        call this%filzC%filter3(zbuff1,zbuff2,this%nxZ, this%nyZ)
        call transpose_z_to_y(zbuff1,this%vhat, this%sp_gpC)

        call transpose_y_to_z(this%what,zbuff3, this%sp_gpE)
        call this%filzE%filter3(zbuff3,zbuff4,this%nxZ, this%nyZ)
        call transpose_z_to_y(zbuff4,this%what, this%sp_gpE)

        if (this%isStratified .or. this%initspinup) then
            call transpose_y_to_z(this%That,zbuff1, this%sp_gpC)
            call this%filzC%filter3(zbuff1,zbuff2,this%nxZ, this%nyZ)
            call transpose_z_to_y(zbuff1,this%That, this%sp_gpC)
        end if

        nullify(zbuff1, zbuff2, zbuff3, zbuff4)
    end subroutine

    subroutine compute_Sijmean(this, Stmp)
        class(igrid), intent(inout), target :: this
        real(rkind),    dimension(this%gpC%xsz(1),   this%gpC%xsz(2),   this%gpC%xsz(3),   6), intent(out), target :: Stmp

        complex(rkind), dimension(this%sp_gpC%ysz(1),this%sp_gpC%ysz(2),this%sp_gpC%ysz(3),3), target :: SfCtmp
        complex(rkind), dimension(:,:,:), pointer :: u_mean3Dhat, v_mean3Dhat, w_mean3Dhat, ctmpz1, ctmpz2, ctmpy1
        real(rkind),    dimension(:,:,:), pointer :: S11_mean3D, S12_mean3D, S13_mean3D, S22_mean3D, S23_mean3D, S33_mean3D, rbuff1
        real(rkind) :: tidSUMreal

        u_mean3Dhat => SfCtmp(:,:,:,1); v_mean3Dhat => SfCtmp(:,:,:,2); w_mean3Dhat => SfCtmp(:,:,:,3)
        S11_mean3D  => Stmp(:,:,:,1);   S12_mean3D  => Stmp(:,:,:,2);   S13_mean3D  => Stmp(:,:,:,3) 
                                        S22_mean3D  => Stmp(:,:,:,4);   S23_mean3D  => Stmp(:,:,:,5) 
                                                                        S33_mean3D  => Stmp(:,:,:,6) 
        ctmpy1 => this%cbuffyC(:,:,:,1); ctmpz1  => this%cbuffzC(:,:,:,1)
        ctmpz2 => this%cbuffzE(:,:,:,1); rbuff1 => this%rbuffxC(:,:,:,1);


        tidSUMreal = real(this%tidSUM, rkind)

        ! compute forward transforms of mean velocities
        call this%spectC%fft(this%u_mean3D/tidSUMreal, u_mean3Dhat)
        call this%spectC%fft(this%v_mean3D/tidSUMreal, v_mean3Dhat)
        call this%spectC%fft(this%w_mean3D/tidSUMreal, w_mean3Dhat)

        ! dudx
        call this%spectC%mTimes_ik1_oop(u_mean3Dhat, ctmpy1)
        call this%spectC%ifft(ctmpy1, S11_mean3D)
         
        ! dudy
        call this%spectC%mTimes_ik2_oop(u_mean3Dhat, ctmpy1)
        call this%spectC%ifft(ctmpy1, S12_mean3D)
         
        ! dvdx
        call this%spectC%mTimes_ik1_oop(v_mean3Dhat, ctmpy1)
        call this%spectC%ifft(ctmpy1, rbuff1)
        S12_mean3D = half*(S12_mean3D + rbuff1)
         
        ! dvdy
        call this%spectC%mTimes_ik2_oop(v_mean3Dhat, ctmpy1)
        call this%spectC%ifft(ctmpy1, S22_mean3D)
         
        ! dwdx
        call this%spectC%mTimes_ik1_oop(w_mean3Dhat, ctmpy1)
        call this%spectC%ifft(ctmpy1, S13_mean3D)
         
        ! dwdy
        call this%spectC%mTimes_ik2_oop(w_mean3Dhat, ctmpy1)
        call this%spectC%ifft(ctmpy1, S23_mean3D)
        
        ! dudz 
        call transpose_y_to_z(u_mean3Dhat, ctmpz1, this%sp_gpC)
        call this%Pade6opZ%ddz_C2E(ctmpz1, ctmpz2, uBC_bottom, uBC_top)
        call this%Pade6opZ%interpz_E2C(ctmpz2, ctmpz1, dUdzBC_bottom, dUdzBC_top)
        call transpose_z_to_y(ctmpz1, u_mean3Dhat, this%sp_gpC)
        call this%spectC%ifft(u_mean3Dhat, rbuff1)
        S13_mean3D = half*(S13_mean3D + rbuff1)

        ! dvdz 
        call transpose_y_to_z(v_mean3Dhat, ctmpz1, this%sp_gpC)
        call this%Pade6opZ%ddz_C2E(ctmpz1, ctmpz2, vBC_bottom, vBC_top)
        call this%Pade6opZ%interpz_E2C(ctmpz2, ctmpz1, dVdzBC_bottom, dVdzBC_top)
        call transpose_z_to_y(ctmpz1, v_mean3Dhat, this%sp_gpC)
        call this%spectC%ifft(v_mean3Dhat, rbuff1)
        S23_mean3D = half*(S23_mean3D + rbuff1)

        ! dwdz 
        call transpose_y_to_z(w_mean3Dhat, ctmpz1, this%sp_gpC)
        call this%Pade6opZ%ddz_C2E(ctmpz1, ctmpz2, wBC_bottom, wBC_top)
        call this%Pade6opZ%interpz_E2C(ctmpz2, ctmpz1, dWdzBC_bottom, dWdzBC_top)
        call transpose_z_to_y(ctmpz1, w_mean3Dhat, this%sp_gpC)
        call this%spectC%ifft(w_mean3Dhat, S33_mean3D)

        nullify(u_mean3Dhat, v_mean3Dhat, w_mean3Dhat, rbuff1, ctmpy1, ctmpz1, ctmpz2, S11_mean3D, S12_mean3D, S13_mean3D, S22_mean3D, S23_mean3D, S33_mean3D)

   end subroutine 


    subroutine compute_vorticity(this)
         class(igrid), intent(inout), target :: this
         real(rkind),    dimension(:,:,:), pointer :: dudx , dudy , dudzC
         real(rkind),    dimension(:,:,:), pointer :: dvdx , dvdy , dvdzC
         real(rkind),    dimension(:,:,:), pointer :: dwdxC, dwdyC, dwdz
        
         
         dudx  => this%duidxjC(:,:,:,1); dudy  => this%duidxjC(:,:,:,2); dudzC => this%duidxjC(:,:,:,3) 
         dvdx  => this%duidxjC(:,:,:,4); dvdy  => this%duidxjC(:,:,:,5); dvdzC => this%duidxjC(:,:,:,6) 
         dwdxC => this%duidxjC(:,:,:,7); dwdyC => this%duidxjC(:,:,:,8); dwdz  => this%duidxjC(:,:,:,9) 


         this%ox = dwdyC - dvdzC
         this%oy = dudzC - dwdxC
         this%oz = dvdx  - dudy

    end subroutine 
    
    subroutine compute_duidxj(this)
        class(igrid), intent(inout), target :: this
        complex(rkind), dimension(:,:,:), pointer :: ctmpz1, ctmpz2
        complex(rkind), dimension(:,:,:), pointer :: ctmpz3, ctmpz4
        complex(rkind), dimension(:,:,:), pointer :: ctmpy1, ctmpy2
        real(rkind),    dimension(:,:,:), pointer :: dudx, dudy, dudz
        real(rkind),    dimension(:,:,:), pointer :: dvdx, dvdy, dvdz
        real(rkind),    dimension(:,:,:), pointer :: dwdx, dwdy, dwdz
        real(rkind),    dimension(:,:,:), pointer :: dvdzC, dudzC 
        real(rkind),    dimension(:,:,:), pointer :: dwdxC, dwdyC
        real(rkind),    dimension(:,:,:), pointer :: dwdzE, dudxE
        real(rkind),    dimension(:,:,:), pointer :: dudyE, dvdxE 
        real(rkind),    dimension(:,:,:), pointer :: dvdyE
        
        complex(rkind), dimension(:,:,:), pointer :: dudxH, dudyH, dudzH 
        complex(rkind), dimension(:,:,:), pointer :: dvdxH, dvdyH, dvdzH
        complex(rkind), dimension(:,:,:), pointer :: dwdxH, dwdyH, dwdzH
        
        complex(rkind), dimension(:,:,:), pointer :: dudxEH, dudyEH, dudzEH 
        complex(rkind), dimension(:,:,:), pointer :: dvdxEH, dvdyEH, dvdzEH
        complex(rkind), dimension(:,:,:), pointer :: dwdxEH, dwdyEH, dwdzEH

        dudx  => this%duidxjC(:,:,:,1); dudy  => this%duidxjC(:,:,:,2); dudzC => this%duidxjC(:,:,:,3) 
        dvdx  => this%duidxjC(:,:,:,4); dvdy  => this%duidxjC(:,:,:,5); dvdzC => this%duidxjC(:,:,:,6) 
        dwdxC => this%duidxjC(:,:,:,7); dwdyC => this%duidxjC(:,:,:,8); dwdz  => this%duidxjC(:,:,:,9) 

        dudxH => this%duidxjChat(:,:,:,1); dudyH => this%duidxjChat(:,:,:,2); dudzH => this%duidxjChat(:,:,:,3) 
        dvdxH => this%duidxjChat(:,:,:,4); dvdyH => this%duidxjChat(:,:,:,5); dvdzH => this%duidxjChat(:,:,:,6) 
        dwdxH => this%duidxjChat(:,:,:,7); dwdyH => this%duidxjChat(:,:,:,8); dwdzH => this%duidxjChat(:,:,:,9) 
        
        dudxEH => this%duidxjEhat(:,:,:,1); dudyEH => this%duidxjEhat(:,:,:,2); dudzEH => this%duidxjEhat(:,:,:,3) 
        dvdxEH => this%duidxjEhat(:,:,:,4); dvdyEH => this%duidxjEhat(:,:,:,5); dvdzEH => this%duidxjEhat(:,:,:,6) 
        dwdxEH => this%duidxjEhat(:,:,:,7); dwdyEH => this%duidxjEhat(:,:,:,8); dwdzEH => this%duidxjEhat(:,:,:,9)
       
        dudxE => this%duidxjE(:,:,:,1); dudyE => this%duidxjE(:,:,:,2); dudz  => this%duidxjE(:,:,:,3)
        dvdxE => this%duidxjE(:,:,:,4); dvdyE => this%duidxjE(:,:,:,5); dvdz  => this%duidxjE(:,:,:,6)
        dwdx  => this%duidxjE(:,:,:,7); dwdy  => this%duidxjE(:,:,:,8); dwdzE => this%duidxjE(:,:,:,9)

        ctmpz1 => this%cbuffzC(:,:,:,1); ctmpz2 => this%cbuffzE(:,:,:,1); 
        ctmpz3 => this%cbuffzC(:,:,:,2); ctmpz4 => this%cbuffzE(:,:,:,2)
        ctmpy1 => this%cbuffyC(:,:,:,1); ctmpy2 => this%cbuffyE(:,:,:,1)

      
        ! dudx
        call this%spectC%mTimes_ik1_oop(this%uhat,dudxH)
        call this%spectC%ifft(dudxH,dudx)
        call this%spectE%mTimes_ik1_oop(this%uEhat,dudxEH)
        call this%spectE%ifft(dudxEH,dudxE)
         
        ! dudy
        call this%spectC%mTimes_ik2_oop(this%uhat,dudyH)
        call this%spectC%ifft(dudyH,dudy)
        call this%spectE%mTimes_ik2_oop(this%uEhat,dudyEH)
        call this%spectE%ifft(dudyEH,dudyE)

        ! dvdx 
        call this%spectC%mTimes_ik1_oop(this%vhat,dvdxH)
        call this%spectC%ifft(dvdxH,dvdx)
        call this%spectE%mTimes_ik1_oop(this%vEhat,dvdxEH)
        call this%spectE%ifft(dvdxEH,dvdxE)

        ! dvdy
        call this%spectC%mTimes_ik2_oop(this%vhat,dvdyH)
        call this%spectC%ifft(dvdyH,dvdy)
        call this%spectE%mTimes_ik2_oop(this%vEhat,dvdyEH)
        call this%spectE%ifft(dvdyEH,dvdyE)

        ! dwdx
        call this%spectC%mTimes_ik1_oop(this%whatC,dwdxH)
        call this%spectC%ifft(dwdxH,dwdxC)
        call this%spectE%mTimes_ik1_oop(this%what, dwdxEH)
        call this%spectE%ifft(dwdxEH,dwdx)

        ! dwdy
        call this%spectC%mTimes_ik2_oop(this%whatC,dwdyH)
        call this%spectC%ifft(dwdyH,dwdyC)
        call this%spectE%mTimes_ik2_oop(this%what,dwdyEH)
        call this%spectE%ifft(dwdyEH,dwdy)
       
        ! dwdz
        call transpose_y_to_z(this%what,ctmpz2,this%sp_gpE)
        call this%Pade6opZ%ddz_E2C(ctmpz2,ctmpz1,wBC_bottom,wBC_top)
        call transpose_z_to_y(ctmpz1,dwdzH,this%sp_gpC)
        call this%spectC%ifft(dwdzH,dwdz)
        call this%Pade6opZ%interpz_C2E(ctmpz1,ctmpz4,dWdzBC_bottom, dWdzBC_top)
        call transpose_z_to_y(ctmpz4,dwdzEH,this%sp_gpE)
        call this%spectE%ifft(dwdzEH,dwdzE)

        ! d2wdz2
        if(.not. this%isinviscid) then
           call this%Pade6opZ%d2dz2_E2E(ctmpz2,ctmpz4,wBC_bottom,wBC_top)
           call transpose_z_to_y(ctmpz4,this%d2wdz2hatE,this%sp_gpE)
        end if

        ! dudz and d2udz2
        call transpose_y_to_z(this%uhat,ctmpz1,this%sp_gpC)
        call this%Pade6opZ%ddz_C2E(ctmpz1,ctmpz2,uBC_bottom,uBC_top)
        call transpose_z_to_y(ctmpz2,dudzEH,this%sp_gpE)
        call this%spectE%ifft(dudzEH,dudz)
        if (.not. this%isinviscid) then
               if ((uBC_top == 0) .or. (uBC_bottom == 0)) then
                  call this%Pade6opZ%ddz_C2E(ctmpz1,ctmpz4,uBC_bottom,uBC_top)
                  call this%Pade6opZ%ddz_E2C(ctmpz4,ctmpz3,dUdzBC_bottom,dUdzBC_top)
               else
                  call this%Pade6opZ%d2dz2_C2C(ctmpz1,ctmpz3,uBC_bottom,uBC_top)
               end if
               call transpose_z_to_y(ctmpz3,this%d2udz2hatC,this%sp_gpC)
        end if 
        call this%Pade6opZ%interpz_E2C(ctmpz2,ctmpz1,dUdzBC_bottom,dUdzBC_top)
        call transpose_z_to_y(ctmpz1,dudzH,this%sp_gpC)
        call this%spectC%ifft(dudzH,dudzC)
      
        ! dvdz and d2vdz2
        call transpose_y_to_z(this%vhat,ctmpz1,this%sp_gpC)
        call this%Pade6opZ%ddz_C2E(ctmpz1,ctmpz2,vBC_bottom,vBC_top)
        call transpose_z_to_y(ctmpz2,dvdzEH,this%sp_gpE)
        call this%spectE%ifft(dvdzEH,dvdz)
        if (.not. this%isinviscid) then
            if ((vBC_top == 0) .or. (vBC_bottom == 0)) then
               call this%Pade6opZ%ddz_C2E(ctmpz1,ctmpz4,vBC_bottom,vBC_top)
               call this%Pade6opZ%ddz_E2C(ctmpz4,ctmpz3,dVdzBC_bottom,dVdzBC_top)
            else
               call this%Pade6opZ%d2dz2_C2C(ctmpz1,ctmpz3,vBC_bottom,vBC_top)
            end if
            call transpose_z_to_y(ctmpz3,this%d2vdz2hatC,this%sp_gpC)
        end if 
        call this%Pade6opZ%interpz_E2C(ctmpz2,ctmpz1,dVdzBC_bottom,dVdzBC_top)
        call transpose_z_to_y(ctmpz1,dvdzH,this%sp_gpC)
        call this%spectC%ifft(dvdzH,dvdzC)

    end subroutine


    subroutine compute_dTdxi(this)
        class(igrid), intent(inout), target :: this
        complex(rkind), dimension(:,:,:), pointer :: ctmpz1, ctmpz2, ctmpz3
        complex(rkind), dimension(:,:,:), pointer :: ctmpy1

        ctmpz1 => this%cbuffzC(:,:,:,1); ctmpz2 => this%cbuffzE(:,:,:,1); 
        ctmpy1 => this%cbuffyC(:,:,:,1); ctmpz3 => this%cbuffzC(:,:,:,2)

        call this%spectC%mtimes_ik1_oop(this%That,this%dTdxH)
        call this%spectC%ifft(this%dTdxH,this%dTdxC)

        call this%spectC%mtimes_ik2_oop(this%That,this%dTdyH)
        call this%spectC%ifft(this%dTdyH,this%dTdyC)
   
        call transpose_y_to_z(this%That, ctmpz1, this%sp_gpC)
        call this%Pade6opZ%ddz_C2E(ctmpz1,ctmpz2,TBC_bottom,TBC_top)
        if (.not. this%isInviscid) then
            call this%Pade6opZ%d2dz2_C2C(ctmpz1,ctmpz3,TBC_bottom, TBC_top)    
            call transpose_z_to_y(ctmpz3,this%d2Tdz2hatC,this%sp_gpC)
        end if 
        call this%Pade6opZ%interpz_E2C(ctmpz2,ctmpz1,dTdzBC_bottom,dTdzBC_top)

        call transpose_z_to_y(ctmpz2, this%dTdzH, this%sp_gpE)
        call this%spectE%ifft(this%dTdzH,this%dTdzE)
        
        call transpose_z_to_y(ctmpz1,this%dTdzHC,this%sp_gpC)
        call this%spectC%ifft(this%dTdzHC,this%dTdzC)



    end subroutine
    
    
    subroutine readRestartFile(this, tid, rid)
        use decomp_2d_io
        use mpi
        use exits, only: message
        use kind_parameters, only: mpirkind
        class(igrid), intent(inout) :: this
        integer, intent(in) :: tid, rid
        character(len=clen) :: tempname, fname
        integer :: ierr, fid 

        write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",rid, "_u.",tid
        fname = this%InputDir(:len_trim(this%InputDir))//"/"//trim(tempname)
        call decomp_2d_read_one(1,this%u,fname, this%gpC)

        write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",rid, "_v.",tid
        fname = this%InputDir(:len_trim(this%InputDir))//"/"//trim(tempname)
        call decomp_2d_read_one(1,this%v,fname, this%gpC)

        write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",rid, "_w.",tid
        fname = this%InputDir(:len_trim(this%InputDir))//"/"//trim(tempname)
        call decomp_2d_read_one(1,this%w,fname, this%gpE)

        if (this%isStratified) then
            write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",rid, "_T.",tid
            fname = this%InputDir(:len_trim(this%InputDir))//"/"//trim(tempname)
            call decomp_2d_read_one(1,this%T,fname, this%gpC)
        end if 

        if (nrank == 0) then
            write(tempname,"(A7,A4,I2.2,A6,I6.6)") "RESTART", "_Run",rid, "_info.",tid
            fname = this%InputDir(:len_trim(this%InputDir))//"/"//trim(tempname)
            fid = 10
            open(unit=fid,file=trim(fname),status="old",action="read")
            read (fid, "(100g15.5)")  this%tsim
            close(fid)
        end if 

        call mpi_barrier(mpi_comm_world, ierr)
        call mpi_bcast(this%tsim,1,mpirkind,0,mpi_comm_world,ierr)
        call mpi_barrier(mpi_comm_world, ierr)
        call message("================= RESTART FILE USED ======================")
        call message(0, "Simulation Time at restart:", this%tsim)
        call message("=================================== ======================")

    end subroutine

    subroutine dumpRestartFile(this)
        use decomp_2d_io
        use mpi
        use exits, only: message
        class(igrid), intent(in) :: this
        character(len=clen) :: tempname, fname
        integer :: ierr

        write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",this%runID, "_u.",this%step
        fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
        call decomp_2d_write_one(1,this%u,fname, this%gpC)

        write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",this%runID, "_v.",this%step
        fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
        call decomp_2d_write_one(1,this%v,fname, this%gpC)

        write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",this%runID, "_w.",this%step
        fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
        call decomp_2d_write_one(1,this%w,fname, this%gpE)

        if (this%isStratified .or. this%initspinup) then
            write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",this%runID, "_T.",this%step
            fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
            call decomp_2d_write_one(1,this%T,fname, this%gpC)
        end if 

        if (nrank == 0) then
            write(tempname,"(A7,A4,I2.2,A6,I6.6)") "RESTART", "_Run",this%runID, "_info.",this%step
            fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
            OPEN(UNIT=10, FILE=trim(fname))
            write(10,"(100g15.5)") this%tsim
            close(10)
        end if 

        call mpi_barrier(mpi_comm_world, ierr)
        call message(1, "Just Dumped a RESTART file")

    end subroutine 


    ! NOTE: If you want to dump an edge field, you need to call in dumpFullField
    ! routine with this%gpE passed in as the 3rd argument. If it's a cell field,
    ! then you don't need to pass in any gp since the default gp is this%gpC
    subroutine dumpFullField(this,arr,label,gp2use)
        use decomp_2d_io
        use mpi
        use exits, only: message
        class(igrid), intent(in) :: this
        character(len=clen) :: tempname, fname
        real(rkind), dimension(:,:,:), intent(in) :: arr
        character(len=4), intent(in) :: label
        type(decomp_info), intent(in), optional :: gp2use

         write(tempname,"(A3,I2.2,A1,A4,A2,I6.6,A4)") "Run",this%runID, "_",label,"_t",this%step,".out"
         fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
         if (present(gp2use)) then
            call decomp_2d_write_one(1,arr,fname,gp2use)
         else
            call decomp_2d_write_one(1,arr,fname,this%gpC)
         end if

    end subroutine


    subroutine dumpVisualizationInfo(this)
        class(igrid), intent(in) :: this
        character(len=clen) :: tempname, fname

        
      if (nrank == 0) then
            write(tempname,"(A3,I2.2,A1,A4,A2,I6.6,A4)") "Run",this%runID, "_","info","_t",this%step,".out"
            fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
            OPEN(UNIT=10, FILE=trim(fname))
            write(10,"(100g17.9)") this%tsim
            write(10,"(100g17.9)") real(this%nx,rkind)
            write(10,"(100g17.9)") real(this%ny,rkind)
            write(10,"(100g17.9)") real(this%nz,rkind)
            close(10)
        end if 
    end subroutine 
    !! STATISTICS !!

    !--------------------------------Beginning 3D Statistics----------------------------------------------
    subroutine init_stats3D(this)
        use exits, only: message
        class(igrid), intent(inout), target :: this
        integer :: nstatsvar, nhorzavgvars, nstv, nhzv


        nstatsvar = 12; nhorzavgvars = 17
        if(this%fastCalcPressure .or. this%storePressure) then
           nstatsvar = nstatsvar + 4
           nhorzavgvars = nhorzavgvars + 3
        endif
        if(.not. this%isInviscid) then
           nstatsvar = nstatsvar + 4
           nhorzavgvars = nhorzavgvars + 4
        endif
        if(this%useSGS) then
           nstatsvar = nstatsvar + 10
           nhorzavgvars = nhorzavgvars + 10
        endif
        if(this%useWindTurbines) then
           nstatsvar = nstatsvar + 4
           nhorzavgvars = nhorzavgvars + 5
        endif
        if(this%isStratified) then
           nstatsvar = nstatsvar + 5
           nhorzavgvars = nhorzavgvars + 5
           if(this%useSGS) then
              nstatsvar = nstatsvar + 3
              nhorzavgvars = nhorzavgvars + 3
           endif
        endif

        allocate(this%debugavg(5),this%debuginst(5))
        allocate(this%inst_horz_avg(5)) ! [ustar, uw, vw, Linv, wT]
        allocate(this%runningSum_sc(5))
        allocate(this%stats3D(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3),nstatsvar))
        allocate(this%horzavgstats(this%nz,nhorzavgvars))

        if(this%useWindTurbines) then
            allocate(this%inst_horz_avg_turb(8*this%WindTurbineArr%nTurbines))
            allocate(this%runningSum_sc_turb(8*this%WindTurbineArr%nTurbines))
            allocate(this%runningSum_turb   (8*this%WindTurbineArr%nTurbines))
        endif

        if (this%computeSpectra) then
            allocate(this%xspectra_mean(this%sp_gpC%ysz(1),this%sp_gpC%ysz(2),this%sp_gpC%ysz(3)))   ! ensure that number of variables for which spectrum is to be computed is smaller than nyg
        end if 

        ! mean velocities
        this%u_mean3D => this%stats3D(:,:,:,1);  this%v_mean3D  => this%stats3D(:,:,:,2);  this%w_mean3D => this%stats3D(:,:,:,3) 
        ! mean squared velocities
        this%uu_mean3D => this%stats3D(:,:,:,4); this%uv_mean3D => this%stats3D(:,:,:,5); this%uw_mean3D => this%stats3D(:,:,:,6)
                                                 this%vv_mean3D => this%stats3D(:,:,:,7); this%vw_mean3D => this%stats3D(:,:,:,8) 
                                                                                          this%ww_mean3D => this%stats3D(:,:,:,9)

        ! triple product of velocities
        this%tketurbtranspx_mean3D => this%stats3D(:,:,:,10)
        this%tketurbtranspy_mean3D => this%stats3D(:,:,:,11)
        this%tketurbtranspz_mean3D => this%stats3D(:,:,:,12)
        nstv = 3 + 6 + 3

        if(this%fastCalcPressure .or. this%storePressure) then
           this%p_mean3D  => this%stats3D(:,:,:,nstv+1)
           this%pu_mean3D => this%stats3D(:,:,:,nstv+2)
           this%pv_mean3D => this%stats3D(:,:,:,nstv+3)
           this%pw_mean3D => this%stats3D(:,:,:,nstv+4)
           nstv = nstv + 4
        endif

        if(.not. this%isInviscid) then
           this%viscdisp_mean3D => this%stats3D(:,:,:,nstv+1)
           this%Siju1_mean3D    => this%stats3D(:,:,:,nstv+2)
           this%Siju2_mean3D    => this%stats3D(:,:,:,nstv+3)
           this%Siju3_mean3D    => this%stats3D(:,:,:,nstv+4)
           nstv = nstv + 4
        endif

        if(this%useSGS)  then
           ! SGS stresses
           this%tau11_mean3D => this%stats3D(:,:,:,nstv+1); this%tau12_mean3D => this%stats3D(:,:,:,nstv+2); this%tau13_mean3D => this%stats3D(:,:,:,nstv+3)
                                                            this%tau22_mean3D => this%stats3D(:,:,:,nstv+4); this%tau23_mean3D => this%stats3D(:,:,:,nstv+5) 
                                                                                                             this%tau33_mean3D => this%stats3D(:,:,:,nstv+6)
           ! SGS dissipation
           this%sgsdissp_mean3D => this%stats3D(:,:,:,nstv+7)

           this%tauu1_mean3D    => this%stats3D(:,:,:,nstv+8)
           this%tauu2_mean3D    => this%stats3D(:,:,:,nstv+9)
           this%tauu3_mean3D    => this%stats3D(:,:,:,nstv+10)

           nstv = nstv + 10

        endif

        if(this%useWindTurbines)  then
           this%turbfx_mean3D => this%stats3D(:,:,:,nstv+1)
           this%turbfy_mean3D => this%stats3D(:,:,:,nstv+2)
           this%turbfz_mean3D => this%stats3D(:,:,:,nstv+3)
           this%uturbf_mean3D => this%stats3D(:,:,:,nstv+4)
           nstv = nstv + 4
        endif

        if (this%isStratified) then
           this%TT_mean3D => this%stats3D(:,:,:,nstv+1);  this%wT_mean3D => this%stats3D(:,:,:,nstv+2);  this%vT_mean3D => this%stats3D(:,:,:,nstv+3)
           this%uT_mean3D => this%stats3D(:,:,:,nstv+4);  this%T_mean3D  => this%stats3D(:,:,:,nstv+5);
           nstv = nstv + 5
           if(this%useSGS) then
              this%q1_mean3D => this%stats3D(:,:,:,nstv+1);  this%q2_mean3D => this%stats3D(:,:,:,nstv+2);  this%q3_mean3D => this%stats3D(:,:,:,nstv+3)
              nstv = nstv + 3
           endif
        end if 

        ! horizontal averages
        ! mean velocities
        this%u_mean   => this%horzavgstats(:,1);  this%v_mean    => this%horzavgstats(:,2);  this%w_mean   => this%horzavgstats(:,3) 
        ! mean squared velocities
        this%uu_mean => this%horzavgstats(:,4); this%uv_mean => this%horzavgstats(:,5); this%uw_mean => this%horzavgstats(:,6)
                                                this%vv_mean => this%horzavgstats(:,7); this%vw_mean => this%horzavgstats(:,8) 
                                                                                        this%ww_mean => this%horzavgstats(:,9)
        ! Dispersive stresses
        this%disperuw_mean => this%horzavgstats(:,10); this%dispervw_mean => this%horzavgstats(:,11)

        ! Mean Kinetic Energy Equation: advection, turbulent transport, dissipation
        this%mkeadv_mean => this%horzavgstats(:,12); this%mkett_mean => this%horzavgstats(:,13); this%mkedisp_mean => this%horzavgstats(:,14)
        ! Turbulent Kinetic Energy Equation: advection, turbulent transport, shear production
        this%tkeadv_mean => this%horzavgstats(:,15); this%tkett_mean => this%horzavgstats(:,16); this%tkeprod_mean => this%horzavgstats(:,17)
        nhzv = 17

        if(this%fastCalcPressure .or. this%storePressure) then
           this%p_mean  => this%horzavgstats(:,nhzv+1)
           this%mkept_mean => this%horzavgstats(:,nhzv+2)
           this%tkept_mean => this%horzavgstats(:,nhzv+3)
           nhzv = nhzv + 3
        endif

        if(.not. this%isInviscid) then
           this%mkevdif_mean => this%horzavgstats(:,nhzv+1)
           this%mkevdsp_mean => this%horzavgstats(:,nhzv+2)
           this%tkevdif_mean => this%horzavgstats(:,nhzv+3)
           this%tkevdsp_mean => this%horzavgstats(:,nhzv+4)
           nhzv = nhzv + 4
        endif

        if(this%useSGS) then
           ! SGS stresses
           this%tau11_mean => this%horzavgstats(:,nhzv+1); this%tau12_mean => this%horzavgstats(:,nhzv+2); this%tau13_mean => this%horzavgstats(:,nhzv+3)
                                                           this%tau22_mean => this%horzavgstats(:,nhzv+4); this%tau23_mean => this%horzavgstats(:,nhzv+5) 
                                                                                                           this%tau33_mean => this%horzavgstats(:,nhzv+6)
           ! SGS dissipation
           this%mkesgst_mean => this%horzavgstats(:,nhzv+7)
           this%mkesgsd_mean => this%horzavgstats(:,nhzv+8)
           this%tkesgst_mean => this%horzavgstats(:,nhzv+9)
           this%tkesgsd_mean => this%horzavgstats(:,nhzv+10)
           nhzv = nhzv + 10

        endif

        if(this%useWindTurbines) then
           this%turbfx_mean => this%horzavgstats(:,nhzv+1)
           this%turbfy_mean => this%horzavgstats(:,nhzv+2)
           this%turbfz_mean => this%horzavgstats(:,nhzv+3)
           this%mketurbf_mean => this%horzavgstats(:,nhzv+4)
           this%tketurbf_mean => this%horzavgstats(:,nhzv+5)
           nhzv = nhzv + 5
        endif

        if (this%isStratified) then
            this%T_mean  => this%horzavgstats(:,nhzv+1);  this%uT_mean => this%horzavgstats(:,nhzv+2);  this%vT_mean => this%horzavgstats(:,nhzv+3)
            this%wT_mean => this%horzavgstats(:,nhzv+4);  this%TT_mean => this%horzavgstats(:,nhzv+5)
            nhzv = nhzv + 5 
           if(this%useSGS) then
              this%q1_mean => this%horzavgstats(:,nhzv+1); this%q2_mean => this%horzavgstats(:,nhzv+2);  this%q3_mean => this%horzavgstats(:,nhzv+3)
              nhzv = nhzv + 3
           endif
        end if 

        this%tidSUM = 0
        this%tprev2 = -1; this%tprev1 = -1;
        this%stats3D = zero
        this%horzavgstats = zero
        this%debugavg = zero
        this%debuginst = zero
        this%inst_horz_avg = zero
        this%runningSum_sc = zero
        if(this%useWindTurbines) then
            this%inst_horz_avg_turb = zero
            this%runningSum_sc_turb = zero
            this%runningSum_turb    = zero
        endif
        this%xspectra_mean = zero

        if((nhzv .ne. nhorzavgvars) .or. (nstv .ne. nstatsvar)) then
            call message(0,"Error in init_stats3D")
            write(*,*) 'nhzv = ', nhzv, nhorzavgvars
            write(*,*) 'nstv = ', nstv, nstatsvar
        endif

        call message(0,"Done init_stats3D")
    end subroutine

    subroutine compute_stats3D(this)
        use kind_parameters, only: mpirkind
        class(igrid), intent(inout), target :: this
        real(rkind), dimension(:,:,:), pointer :: rbuff0, rbuff1, rbuff2, rbuff2E, rbuff3E, rbuff3, rbuff1E, rbuff4
        integer :: j, k, jindx, ierr
        real(rkind),    dimension(:,:,:), pointer :: dudx, dudy
        real(rkind),    dimension(:,:,:), pointer :: dvdx, dvdy
        real(rkind),    dimension(:,:,:), pointer :: dwdz
        real(rkind),    dimension(:,:,:), pointer :: dvdzC, dudzC
        real(rkind),    dimension(:,:,:), pointer :: dwdxC, dwdyC

        rbuff0  => this%rbuffxC(:,:,:,1); rbuff1  => this%rbuffxC(:,:,:,2);
        rbuff2  => this%rbuffyC(:,:,:,1);
        rbuff2E => this%rbuffyE(:,:,:,1); rbuff3E => this%rbuffzE(:,:,:,1);
        rbuff3 => this%rbuffzC(:,:,:,1);  rbuff1E => this%rbuffxE(:,:,:,1)
        rbuff4 => this%rbuffzC(:,:,:,2);

        dudx  => this%duidxjC(:,:,:,1); dudy  => this%duidxjC(:,:,:,2); dudzC => this%duidxjC(:,:,:,3); 
        dvdx  => this%duidxjC(:,:,:,4); dvdy  => this%duidxjC(:,:,:,5); dvdzC => this%duidxjC(:,:,:,6); 
        dwdxC => this%duidxjC(:,:,:,7); dwdyC => this%duidxjC(:,:,:,8); dwdz  => this%duidxjC(:,:,:,9); 

        this%tidSUM = this%tidSUM + 1

        ! compute u*w on E, interpolate to C
        rbuff1E = this%uE * this%w
        call transpose_x_to_y(rbuff1E,rbuff2E,this%gpE)
        call transpose_y_to_z(rbuff2E,rbuff3E,this%gpE)
        call this%OpsPP%InterpZ_Edge2Cell(rbuff3E,rbuff3)
        call transpose_z_to_y(rbuff3,rbuff2,this%gpC)
        call transpose_y_to_x(rbuff2,rbuff1,this%gpC)

        ! compute v*w on E, interpolate to C
        rbuff1E = this%vE * this%w
        call transpose_x_to_y(rbuff1E,rbuff2E,this%gpE)
        call transpose_y_to_z(rbuff2E,rbuff3E,this%gpE)
        call this%OpsPP%InterpZ_Edge2Cell(rbuff3E,rbuff3)
        call transpose_z_to_y(rbuff3,rbuff2,this%gpC)
        call transpose_y_to_x(rbuff2,rbuff0,this%gpC)

        ! Compute u,v,wC - mean
        if(this%normByUstar) then
            this%u_mean3D = this%u_mean3D + this%u/this%sgsmodel%get_ustar()
            this%v_mean3D = this%v_mean3D + this%v/this%sgsmodel%get_ustar()
            this%w_mean3D = this%w_mean3D + this%wC/this%sgsmodel%get_ustar()

            this%uu_mean3D = this%uu_mean3D + this%u * this%u /this%sgsmodel%get_ustar()**2
            this%uv_mean3D = this%uv_mean3D + this%u * this%v /this%sgsmodel%get_ustar()**2
            this%uw_mean3D = this%uw_mean3D + rbuff1          /this%sgsmodel%get_ustar()**2
            this%vv_mean3D = this%vv_mean3D + this%v * this%v /this%sgsmodel%get_ustar()**2
            this%vw_mean3D = this%vw_mean3D + rbuff0          /this%sgsmodel%get_ustar()**2
            this%ww_mean3D = this%ww_mean3D + this%wC* this%wC/this%sgsmodel%get_ustar()**2
        else
            this%u_mean3D = this%u_mean3D + this%u
            this%v_mean3D = this%v_mean3D + this%v
            this%w_mean3D = this%w_mean3D + this%wC

            this%uu_mean3D = this%uu_mean3D + this%u * this%u
            this%uv_mean3D = this%uv_mean3D + this%u * this%v
            this%uw_mean3D = this%uw_mean3D + rbuff1
            this%vv_mean3D = this%vv_mean3D + this%v * this%v
            this%vw_mean3D = this%vw_mean3D + rbuff0
            this%ww_mean3D = this%ww_mean3D + this%wC* this%wC
        endif

        ! triple correlation for transport term in TKE budget -- triple product should be dealiased
        rbuff0 = this%u*this%u + this%v*this%v + this%wC*this%wC
        this%tketurbtranspx_mean3D = this%tketurbtranspx_mean3D + this%u *rbuff0
        this%tketurbtranspy_mean3D = this%tketurbtranspy_mean3D + this%v *rbuff0
        this%tketurbtranspz_mean3D = this%tketurbtranspz_mean3D + this%wC*rbuff0

        if(this%fastCalcPressure .or. this%storePressure) then
            if(this%normByUstar) then
                this%p_mean3D  = this%p_mean3D  + this%pressure          /this%sgsmodel%get_ustar()**2
                this%pu_mean3D = this%pu_mean3D + this%pressure * this%u /this%sgsmodel%get_ustar()**3
                this%pv_mean3D = this%pv_mean3D + this%pressure * this%v /this%sgsmodel%get_ustar()**3
                this%pw_mean3D = this%pw_mean3D + this%pressure * this%wC/this%sgsmodel%get_ustar()**3
            else
                this%p_mean3D  = this%p_mean3D  + this%pressure
                this%pu_mean3D = this%pu_mean3D + this%pressure * this%u
                this%pv_mean3D = this%pv_mean3D + this%pressure * this%v
                this%pw_mean3D = this%pw_mean3D + this%pressure * this%wC
            endif
        endif

        if(.not. this%isInviscid) then
            ! for viscous dissipation
            rbuff0 = dudx*dudx + dvdy*dvdy + dwdz*dwdz &
                   + half*( (dudy + dvdx)**2 + (dudzC + dwdxC)**2 + (dvdzC + dwdyC)**2) ! half here is two/four

            if(this%normByUstar) then
                ! viscous dissipation
                this%viscdisp_mean3D = this%viscdisp_mean3D + rbuff0/(this%sgsmodel%get_ustar()**2)

                ! viscous diffusion
                this%Siju1_mean3D = this%Siju1_mean3D + (dudx*this%u + &
                                                         half*(dudy  + dvdx )*this%v + &
                                                         half*(dudzC + dwdxC)*this%wC)/(this%sgsmodel%get_ustar()**2)

                this%Siju2_mean3D = this%Siju2_mean3D + (half*(dudy  + dvdx )*this%u + &
                                                         dvdy*this%v + &
                                                         half*(dvdzC + dwdyC)*this%wC)/(this%sgsmodel%get_ustar()**2)

                this%Siju3_mean3D = this%Siju3_mean3D + (half*(dudzC + dwdxC)*this%u + &
                                                         half*(dvdzC + dwdyC)*this%v + &
                                                         dwdz*this%wC)/(this%sgsmodel%get_ustar()**2)
            else
                ! viscous dissipation
                this%viscdisp_mean3D = this%viscdisp_mean3D + rbuff0

                ! viscous diffusion
                this%Siju1_mean3D = this%Siju1_mean3D + (dudx*this%u + &
                                                         half*(dudy  + dvdx )*this%v + &
                                                         half*(dudzC + dwdxC)*this%wC)

                this%Siju2_mean3D = this%Siju2_mean3D + (half*(dudy  + dvdx )*this%u + &
                                                         dvdy*this%v + &
                                                         half*(dvdzC + dwdyC)*this%wC)

                this%Siju3_mean3D = this%Siju3_mean3D + (half*(dudzC + dwdxC)*this%u + &
                                                         half*(dvdzC + dwdyC)*this%v + &
                                                         dwdz*this%wC)
            endif
        endif

        if(this%useSGS) then

            !write(300+nrank,'(i4,3(e19.12,1x))') 1, this%inst_horz_avg(2:3)
            call mpi_bcast(this%inst_horz_avg(2:3),2,mpirkind,0,mpi_comm_world,ierr)
            !write(300+nrank,'(i4,3(e19.12,1x))') 2, this%inst_horz_avg(2:3)


            ! interpolate tau13 from E to C
            call transpose_x_to_y(this%tau13,rbuff2E,this%gpE)
            call transpose_y_to_z(rbuff2E,rbuff3E,this%gpE) 
            !if(nrank==0) then
            !    write(*,*) '-----------------'
            !    write(*,*) rbuff3E(1,1,1:2)
            !endif
            !write(200+nrank,*) 1, this%tsim, this%inst_horz_avg(2)
            !rbuff3E(:,:,1) = this%sgsmodel%tauijWMhat_in_Z(:,:,1,1)    !this%inst_horz_avg(2)     !--- =-(this%sgsmodel%get_ustar()**2) is correct only for Moeng's Wall Model, not for Bou-Zeid's model
            !if(nrank==0) then
            !    write(*,*) rbuff3E(1,1,1:2)
            !endif
            !this%debuginst(1) = p_sum(sum(rbuff3E(:,:,1)))*this%meanFact
            !write(200+nrank,*) 2, this%tsim, this%debuginst(1)
            !this%debuginst(2) = p_sum(sum(rbuff3E(:,:,2)))*this%meanFact
            !this%debuginst(3) = p_sum(sum(rbuff3E(:,:,3)))*this%meanFact
            call this%OpsPP%InterpZ_Edge2Cell(rbuff3E,rbuff3)
            !this%debuginst(4) = p_sum(sum(rbuff3(:,:,1)))*this%meanFact
            !this%debuginst(5) = p_sum(sum(rbuff3(:,:,2)))*this%meanFact
            !this%debugavg(:) = this%debugavg(:) + this%debuginst(:)
            !if(nrank==0) then
            !    write(nrank+100,'(11(e19.12,1x))') this%tsim, this%debuginst, this%debugavg/real(this%tidSUM, rkind)
            !    !write(*,*) '-----------------'
            !endif
            call transpose_z_to_y(rbuff3,rbuff2,this%gpC)
            call transpose_y_to_x(rbuff2,this%tauSGS_ij(:,:,:,3),this%gpC)

            ! interpolate tau23 from E to C
            call transpose_x_to_y(this%tau23,rbuff2E,this%gpE)
            call transpose_y_to_z(rbuff2E,rbuff3E,this%gpE)
            !rbuff3E(:,:,1) = this%sgsmodel%tauijWMhat_in_Z(:,:,1,2)    !this%inst_horz_avg(3)     !--- =-(this%sgsmodel%get_ustar()**2) is correct only for Moeng's Wall Model, not for Bou-Zeid's model
            call this%OpsPP%InterpZ_Edge2Cell(rbuff3E,rbuff3)
            call transpose_z_to_y(rbuff3,rbuff2,this%gpC)
            call transpose_y_to_x(rbuff2,this%tauSGS_ij(:,:,:,5),this%gpC)

            ! sgs dissipation
            !rbuff1 = this%tauSGS_ij(:,:,:,1)*this%tauSGS_ij(:,:,:,1) + &
            !         this%tauSGS_ij(:,:,:,4)*this%tauSGS_ij(:,:,:,4) + &
            !         this%tauSGS_ij(:,:,:,6)*this%tauSGS_ij(:,:,:,6)
            !rbuff1 = rbuff1 + two*(this%tauSGS_ij(:,:,:,2)*this%tauSGS_ij(:,:,:,2) + &
            !                       this%tauSGS_ij(:,:,:,3)*this%tauSGS_ij(:,:,:,3) + &
            !                       this%tauSGS_ij(:,:,:,5)*this%tauSGS_ij(:,:,:,5) )
            !rbuff1 = rbuff1/(this%nu_SGS + 1.0d-14)         ! note: factor of half is in dump_stats
            rbuff1 = this%tauSGS_ij(:,:,:,1)*dudx + this%tauSGS_ij(:,:,:,4)*dvdy + this%tauSGS_ij(:,:,:,6)*dwdz &
                   + this%tauSGS_ij(:,:,:,2)*(dudy + dvdx) + this%tauSGS_ij(:,:,:,3)*(dudzC + dwdxC) &
                   + this%tauSGS_ij(:,:,:,5)*(dvdzC + dwdyC)

            if(this%normByUstar) then
                this%tau11_mean3D = this%tau11_mean3D + this%tauSGS_ij(:,:,:,1)/(this%sgsmodel%get_ustar()**2)
                this%tau12_mean3D = this%tau12_mean3D + this%tauSGS_ij(:,:,:,2)/(this%sgsmodel%get_ustar()**2)
                this%tau13_mean3D = this%tau13_mean3D + this%tauSGS_ij(:,:,:,3)/(this%sgsmodel%get_ustar()**2)
                this%tau22_mean3D = this%tau22_mean3D + this%tauSGS_ij(:,:,:,4)/(this%sgsmodel%get_ustar()**2)
                this%tau23_mean3D = this%tau23_mean3D + this%tauSGS_ij(:,:,:,5)/(this%sgsmodel%get_ustar()**2)
                this%tau33_mean3D = this%tau33_mean3D + this%tauSGS_ij(:,:,:,6)/(this%sgsmodel%get_ustar()**2)

                ! factor of H in normalization is missing from all statistics below
                this%sgsdissp_mean3D = this%sgsdissp_mean3D + rbuff1/(this%sgsmodel%get_ustar()**3)

                this%tauu1_mean3D = this%tauu1_mean3D + (this%tauSGS_ij(:,:,:,1)*this%u + &
                                                         this%tauSGS_ij(:,:,:,2)*this%v + &
                                                         this%tauSGS_ij(:,:,:,3)*this%wC)/(this%sgsmodel%get_ustar()**3)

                this%tauu2_mean3D = this%tauu2_mean3D + (this%tauSGS_ij(:,:,:,2)*this%u + &
                                                         this%tauSGS_ij(:,:,:,4)*this%v + &
                                                         this%tauSGS_ij(:,:,:,5)*this%wC)/(this%sgsmodel%get_ustar()**3)

                this%tauu3_mean3D = this%tauu3_mean3D + (this%tauSGS_ij(:,:,:,3)*this%u + &
                                                         this%tauSGS_ij(:,:,:,5)*this%v + &
                                                         this%tauSGS_ij(:,:,:,6)*this%wC)/(this%sgsmodel%get_ustar()**3)
            else
                this%tau11_mean3D = this%tau11_mean3D + this%tauSGS_ij(:,:,:,1)
                this%tau12_mean3D = this%tau12_mean3D + this%tauSGS_ij(:,:,:,2)
                this%tau13_mean3D = this%tau13_mean3D + this%tauSGS_ij(:,:,:,3)
                this%tau22_mean3D = this%tau22_mean3D + this%tauSGS_ij(:,:,:,4)
                this%tau23_mean3D = this%tau23_mean3D + this%tauSGS_ij(:,:,:,5)
                this%tau33_mean3D = this%tau33_mean3D + this%tauSGS_ij(:,:,:,6)

                this%sgsdissp_mean3D = this%sgsdissp_mean3D + rbuff1

                this%tauu1_mean3D = this%tauu1_mean3D + (this%tauSGS_ij(:,:,:,1)*this%u + &
                                                         this%tauSGS_ij(:,:,:,2)*this%v + &
                                                         this%tauSGS_ij(:,:,:,3)*this%wC)

                this%tauu2_mean3D = this%tauu2_mean3D + (this%tauSGS_ij(:,:,:,2)*this%u + &
                                                         this%tauSGS_ij(:,:,:,4)*this%v + &
                                                         this%tauSGS_ij(:,:,:,5)*this%wC)

                this%tauu3_mean3D = this%tauu3_mean3D + (this%tauSGS_ij(:,:,:,3)*this%u + &
                                                         this%tauSGS_ij(:,:,:,5)*this%v + &
                                                         this%tauSGS_ij(:,:,:,6)*this%wC)
            endif

        endif

        if(this%useWindTurbines) then
           if(this%normByUstar) then
              this%turbfx_mean3D = this%turbfx_mean3D + this%WindTurbineArr%fx/(this%sgsmodel%get_ustar()**3)
              this%turbfy_mean3D = this%turbfy_mean3D + this%WindTurbineArr%fy/(this%sgsmodel%get_ustar()**3)
              this%turbfz_mean3D = this%turbfz_mean3D + this%WindTurbineArr%fz/(this%sgsmodel%get_ustar()**3)
              this%uturbf_mean3D = this%uturbf_mean3D + (this%u *this%WindTurbineArr%fx + &
                                                         this%v *this%WindTurbineArr%fy + &
                                                         this%wC*this%WindTurbineArr%fz)/(this%sgsmodel%get_ustar()**3)
           else
              this%turbfx_mean3D = this%turbfx_mean3D + this%WindTurbineArr%fx
              this%turbfy_mean3D = this%turbfy_mean3D + this%WindTurbineArr%fy
              this%turbfz_mean3D = this%turbfz_mean3D + this%WindTurbineArr%fz
              this%uturbf_mean3D = this%uturbf_mean3D + this%u *this%WindTurbineArr%fx + &
                                                        this%v *this%WindTurbineArr%fy + &
                                                        this%wC*this%WindTurbineArr%fz
           endif
        endif

        if(this%isStratified) then
            ! compute T*w on E, interpolate to C
            rbuff1E = this%TE * this%w
            call transpose_x_to_y(rbuff1E,rbuff2E,this%gpE)
            call transpose_y_to_z(rbuff2E,rbuff3E,this%gpE)
            call this%OpsPP%InterpZ_Edge2Cell(rbuff3E,rbuff3)
            call transpose_z_to_y(rbuff3,rbuff2,this%gpC)
            call transpose_y_to_x(rbuff2,rbuff0,this%gpC)

            if(this%normByUstar) then
                this%T_mean3D = this%T_mean3D + this%T*this%sgsmodel%get_ustar()/this%wTh_surf
                this%uT_mean3D = this%uT_mean3D + this%T * this%u /this%wTh_surf
                this%vT_mean3D = this%vT_mean3D + this%T * this%v /this%wTh_surf
                this%wT_mean3D = this%wT_mean3D + rbuff0          /this%wTh_surf
                this%TT_mean3D = this%TT_mean3D + this%T * this%T*(this%sgsmodel%get_ustar()/this%wTh_surf)**2
            else
                this%T_mean3D = this%T_mean3D + this%T
                this%uT_mean3D = this%uT_mean3D + this%T * this%u
                this%vT_mean3D = this%vT_mean3D + this%T * this%v
                this%wT_mean3D = this%wT_mean3D + rbuff0
                this%TT_mean3D = this%TT_mean3D + this%T * this%T
            endif
            if(this%useSGS) then
                ! interpolate q3 from E to C
                call transpose_x_to_y(this%q3,rbuff2E,this%gpE)
                call transpose_y_to_z(rbuff2E,rbuff3E,this%gpE)
                rbuff3E(:,:,1) = this%wTh_surf
                call this%OpsPP%InterpZ_Edge2Cell(rbuff3E,rbuff3)
                call transpose_z_to_y(rbuff3,rbuff2,this%gpC)
                call transpose_y_to_x(rbuff2,rbuff1,this%gpC)

                if(this%normByUstar) then
                    this%q1_mean3D = this%q1_mean3D + this%q1/this%wTh_surf
                    this%q2_mean3D = this%q2_mean3D + this%q2/this%wTh_surf
                    this%q3_mean3D = this%q3_mean3D + rbuff1/this%wTh_surf
                else
                    this%q1_mean3D = this%q1_mean3D + this%q1
                    this%q2_mean3D = this%q2_mean3D + this%q2
                    this%q3_mean3D = this%q3_mean3D + rbuff1
                endif
            endif

        endif

        if (this%computeSpectra) then
            ! compute 1D spectra ---- make sure that number of variables for which spectra are computed is smaller than nyg
            ! For each variable, at each y, z, location, x-spectrum is computed first, and then averaged over time and y-direction
            jindx = 1    ! u
            call this%spectC%fft1_x2y(this%u,this%cbuffyC(:,:,:,1))
            do k = 1, size(this%cbuffyC, 3)
              do j = 1, size(this%cbuffyC, 2)
                this%xspectra_mean(:,jindx,k) = this%xspectra_mean(:,jindx,k) + abs(this%cbuffyC(:,j,k,1))
              end do
            end do

            jindx = 2    ! v
            call this%spectC%fft1_x2y(this%v,this%cbuffyC(:,:,:,1))
            do k = 1, size(this%cbuffyC, 3)
              do j = 1, size(this%cbuffyC, 2)
                this%xspectra_mean(:,jindx,k) = this%xspectra_mean(:,jindx,k) + abs(this%cbuffyC(:,j,k,1))
              end do
            end do
            
            jindx = 3    ! w
            call this%spectC%fft1_x2y(this%wC,this%cbuffyC(:,:,:,1))
            do k = 1, size(this%cbuffyC, 3)
              do j = 1, size(this%cbuffyC, 2)
                this%xspectra_mean(:,jindx,k) = this%xspectra_mean(:,jindx,k) + abs(this%cbuffyC(:,j,k,1))
              end do
            end do
            
            jindx = 4    ! KE
            call this%spectC%fft1_x2y(half*(this%u**2+this%v**2+this%wC**2),this%cbuffyC(:,:,:,1))
            do k = 1, size(this%cbuffyC, 3)
              do j = 1, size(this%cbuffyC, 2)
                this%xspectra_mean(:,jindx,k) = this%xspectra_mean(:,jindx,k) + abs(this%cbuffyC(:,j,k,1))
              end do
            end do
            
            if(this%isStratified) then
                jindx = jindx + 1    ! T
                call this%spectC%fft1_x2y(this%T,this%cbuffyC(:,:,:,1))
                do k = 1, size(this%cbuffyC, 3)
                  do j = 1, size(this%cbuffyC, 2)
                    this%xspectra_mean(:,jindx,k) = this%xspectra_mean(:,jindx,k) + abs(this%cbuffyC(:,j,k,1))
                  end do
                end do
            endif
            
            if(this%fastCalcPressure .or. this%storePressure) then
                jindx = jindx + 1    ! p
                call this%spectC%fft1_x2y(this%pressure,this%cbuffyC(:,:,:,1))
                do k = 1, size(this%cbuffyC, 3)
                  do j = 1, size(this%cbuffyC, 2)
                    this%xspectra_mean(:,jindx,k) = this%xspectra_mean(:,jindx,k) + abs(this%cbuffyC(:,j,k,1))
                  end do
                end do
            endif
        end if 

        ! instantaneous horizontal averages of some quantities
        this%inst_horz_avg(1) = this%sgsmodel%get_ustar()
        ! this%inst_horz(2) and (3) are computed on nrank==0 proc in getRHS_SGS_WallM
        ! broadcast to all other procs above in this subroutine
        ! do nothing about inst_horz_avg(2:3) herre
        if(this%isStratified) then
            this%inst_horz_avg(4) = this%invObLength
            this%inst_horz_avg(5) = this%wTh_surf
        endif
        ! this%inst_horz_avg_turb(1:5*this%WindTurbineArr%nTurbines) is computed in this%WindTurbineArr%getForceRHS
        this%runningSum_sc = this%runningSum_sc + this%inst_horz_avg
            !write(200+nrank,*) 3, this%tsim, this%runningSum_sc(2)
        if(this%useWindTurbines) this%runningSum_sc_turb = this%runningSum_sc_turb + this%inst_horz_avg_turb

        nullify(rbuff0,rbuff1,rbuff2,rbuff3,rbuff2E,rbuff3E,rbuff4,rbuff1E)
        nullify(dudx, dudy, dudzC, dvdx, dvdy, dvdzC, dwdxC, dwdyC, dwdz)

    end subroutine

    subroutine Delete_file_if_present(this, tempname)
        class(igrid), intent(inout) :: this
        character(len=clen), intent(in) :: tempname
        character(len=clen) :: fname

        fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
        open(unit=11, file=fname, status='old',iostat=ierr); if(ierr == 0) close(11, status='delete')

    end subroutine

    subroutine DeletePrevStats3DFiles(this)
        class(igrid), intent(inout) :: this
        character(len=clen) :: tempname

        if(nrank==0) then
          ! delete stats files corresponding to tprev2
          write(tempname,"(A3,I2.2,A5,I6.6,A6)") "Run",this%runID, "_um_t",this%tprev2,".3Dstt";  call this%Delete_file_if_present(tempname)
          write(tempname,"(A3,I2.2,A5,I6.6,A6)") "Run",this%runID, "_vm_t",this%tprev2,".3Dstt";  call this%Delete_file_if_present(tempname)
          write(tempname,"(A3,I2.2,A5,I6.6,A6)") "Run",this%runID, "_wm_t",this%tprev2,".3Dstt";  call this%Delete_file_if_present(tempname)
          write(tempname,"(A3,I2.2,A6,I6.6,A6)") "Run",this%runID, "_uum_t",this%tprev2,".3Dstt"; call this%Delete_file_if_present(tempname)
          write(tempname,"(A3,I2.2,A6,I6.6,A6)") "Run",this%runID, "_uvm_t",this%tprev2,".3Dstt"; call this%Delete_file_if_present(tempname)
          write(tempname,"(A3,I2.2,A6,I6.6,A6)") "Run",this%runID, "_uwm_t",this%tprev2,".3Dstt"; call this%Delete_file_if_present(tempname)
          write(tempname,"(A3,I2.2,A6,I6.6,A6)") "Run",this%runID, "_vvm_t",this%tprev2,".3Dstt"; call this%Delete_file_if_present(tempname)
          write(tempname,"(A3,I2.2,A6,I6.6,A6)") "Run",this%runID, "_vwm_t",this%tprev2,".3Dstt"; call this%Delete_file_if_present(tempname)
          write(tempname,"(A3,I2.2,A6,I6.6,A6)") "Run",this%runID, "_wwm_t",this%tprev2,".3Dstt"; call this%Delete_file_if_present(tempname)
          write(tempname,"(A3,I2.2,A7,I6.6,A6)") "Run",this%runID, "_tktt_t",this%tprev2,".3Dstt"; call this%Delete_file_if_present(tempname)
          write(tempname,"(A3,I2.2,A7,I6.6,A6)") "Run",this%runID, "_tkma_t",this%tprev2,".3Dstt"; call this%Delete_file_if_present(tempname)
          write(tempname,"(A3,I2.2,A7,I6.6,A6)") "Run",this%runID, "_tkpr_t",this%tprev2,".3Dstt"; call this%Delete_file_if_present(tempname)
          write(tempname,"(A3,I2.2,A7,I6.6,A6)") "Run",this%runID, "_mktt_t",this%tprev2,".3Dstt"; call this%Delete_file_if_present(tempname)
          write(tempname,"(A3,I2.2,A7,I6.6,A6)") "Run",this%runID, "_mkma_t",this%tprev2,".3Dstt"; call this%Delete_file_if_present(tempname)
          if(this%fastCalcPressure .or. this%storePressure) then
              write(tempname,"(A3,I2.2,A5,I6.6,A6)") "Run",this%runID, "_pm_t",this%tprev2,".3Dstt";   call this%Delete_file_if_present(tempname)
              write(tempname,"(A3,I2.2,A7,I6.6,A6)") "Run",this%runID, "_mkpt_t",this%tprev2,".3Dstt"; call this%Delete_file_if_present(tempname)
              write(tempname,"(A3,I2.2,A7,I6.6,A6)") "Run",this%runID, "_tkpt_t",this%tprev2,".3Dstt"; call this%Delete_file_if_present(tempname)
          endif
          if(this%useSGS) then
              write(tempname,"(A3,I2.2,A6,I6.6,A6)") "Run",this%runID, "_t11_t",this%tprev2,".3Dstt";    call this%Delete_file_if_present(tempname)
              write(tempname,"(A3,I2.2,A6,I6.6,A6)") "Run",this%runID, "_t12_t",this%tprev2,".3Dstt";    call this%Delete_file_if_present(tempname)
              write(tempname,"(A3,I2.2,A6,I6.6,A6)") "Run",this%runID, "_t13_t",this%tprev2,".3Dstt";    call this%Delete_file_if_present(tempname)
              write(tempname,"(A3,I2.2,A6,I6.6,A6)") "Run",this%runID, "_t22_t",this%tprev2,".3Dstt";    call this%Delete_file_if_present(tempname)
              write(tempname,"(A3,I2.2,A6,I6.6,A6)") "Run",this%runID, "_t23_t",this%tprev2,".3Dstt";    call this%Delete_file_if_present(tempname)
              write(tempname,"(A3,I2.2,A6,I6.6,A6)") "Run",this%runID, "_t33_t",this%tprev2,".3Dstt";    call this%Delete_file_if_present(tempname)
              write(tempname,"(A3,I2.2,A9,I6.6,A6)") "Run",this%runID, "_mksgsd_t",this%tprev2,".3Dstt"; call this%Delete_file_if_present(tempname)
              write(tempname,"(A3,I2.2,A9,I6.6,A6)") "Run",this%runID, "_tksgsd_t",this%tprev2,".3Dstt"; call this%Delete_file_if_present(tempname)
              write(tempname,"(A3,I2.2,A9,I6.6,A6)") "Run",this%runID, "_mksgst_t",this%tprev2,".3Dstt"; call this%Delete_file_if_present(tempname)
              write(tempname,"(A3,I2.2,A9,I6.6,A6)") "Run",this%runID, "_tksgst_t",this%tprev2,".3Dstt"; call this%Delete_file_if_present(tempname)
          endif

          if(this%useWindTurbines) then
              write(tempname,"(A3,I2.2,A7,I6.6,A6)") "Run",this%runID, "_trbx_t",this%tprev2,".3Dstt"; call this%Delete_file_if_present(tempname)
              write(tempname,"(A3,I2.2,A7,I6.6,A6)") "Run",this%runID, "_trby_t",this%tprev2,".3Dstt"; call this%Delete_file_if_present(tempname)
              write(tempname,"(A3,I2.2,A7,I6.6,A6)") "Run",this%runID, "_trbz_t",this%tprev2,".3Dstt"; call this%Delete_file_if_present(tempname)
              write(tempname,"(A3,I2.2,A9,I6.6,A6)") "Run",this%runID, "_mktrbf_t",this%tprev2,".3Dstt"; call this%Delete_file_if_present(tempname)
              write(tempname,"(A3,I2.2,A9,I6.6,A6)") "Run",this%runID, "_tktrbf_t",this%tprev2,".3Dstt"; call this%Delete_file_if_present(tempname)
          endif

          if(this%isStratified) then
              write(tempname,"(A3,I2.2,A5,I6.6,A6)") "Run",this%runID, "_Tm_t",this%tprev2,".3Dstt";  call this%Delete_file_if_present(tempname)
              write(tempname,"(A3,I2.2,A6,I6.6,A6)") "Run",this%runID, "_uTm_t",this%tprev2,".3Dstt"; call this%Delete_file_if_present(tempname)
              write(tempname,"(A3,I2.2,A6,I6.6,A6)") "Run",this%runID, "_vTm_t",this%tprev2,".3Dstt"; call this%Delete_file_if_present(tempname)
              write(tempname,"(A3,I2.2,A6,I6.6,A6)") "Run",this%runID, "_wTm_t",this%tprev2,".3Dstt"; call this%Delete_file_if_present(tempname)
              write(tempname,"(A3,I2.2,A6,I6.6,A6)") "Run",this%runID, "_TTm_t",this%tprev2,".3Dstt"; call this%Delete_file_if_present(tempname)
              if(this%useSGS) then
                write(tempname,"(A3,I2.2,A5,I6.6,A6)") "Run",this%runID, "_q1_t",this%tprev2,".3Dstt"; call this%Delete_file_if_present(tempname)
                write(tempname,"(A3,I2.2,A5,I6.6,A6)") "Run",this%runID, "_q2_t",this%tprev2,".3Dstt"; call this%Delete_file_if_present(tempname)
                write(tempname,"(A3,I2.2,A5,I6.6,A6)") "Run",this%runID, "_q3_t",this%tprev2,".3Dstt"; call this%Delete_file_if_present(tempname)
              endif
          endif
        endif
    end subroutine

    subroutine dump_stats3D(this)
        use basic_io, only: write_2d_ascii
        use decomp_2d_io
        use kind_parameters, only: mpirkind
        class(igrid), intent(inout), target :: this
      ! compute horizontal averages and dump .stt files
      ! overwrite previously written out 3D stats dump
        real(rkind), dimension(:,:,:), pointer :: rbuff0, rbuff1, rbuff2, rbuff3, rbuff4, rbuff5, rbuff6, rbuff3E, rbuff4E
        real(rkind), dimension(:,:,:), pointer :: S11_mean3D, S12_mean3D, S13_mean3D, S22_mean3D, S23_mean3D, S33_mean3D
        complex(rkind), dimension(:,:,:), pointer :: cbuffy1, cbuffy2
        real(rkind), dimension(this%sp_gpC%ysz(1),this%sp_gpC%ysz(2),this%sp_gpC%ysz(3)) :: tmpvar
        real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3),6), target :: Stmp
        character(len=clen) :: tempname, fname
        real(rkind) :: tidSUMreal, normfac
        integer :: tid, dirid, decompdir, jindx, nspectra

        tid = this%step

        ! Ensure only two sets of 3Dstats files are kept
        if(this%tprev2 > 0) then
          call this%DeletePrevStats3DFiles()
        endif
        ! Now update tprev2 and tprev1 
        this%tprev2 = this%tprev1
        this%tprev1 = this%step

        rbuff0 => this%rbuffxC(:,:,:,2);
        rbuff1 => this%rbuffxC(:,:,:,1);  rbuff2 => this%rbuffyC(:,:,:,1)
        rbuff3 => this%rbuffzC(:,:,:,1);  rbuff4 => this%rbuffzC(:,:,:,2)
        rbuff5 => this%rbuffzC(:,:,:,3);  rbuff6 => this%rbuffzC(:,:,:,4)

        cbuffy1 => this%cbuffyC(:,:,:,1); cbuffy2 => this%cbuffyC(:,:,:,2)
        rbuff3E => this%rbuffzE(:,:,:,1); rbuff4E => this%rbuffzE(:,:,:,2);

        tidSUMreal = real(this%tidSUM, rkind)

        ! u_avg
        call transpose_x_to_y(this%u_mean3D/tidSUMreal, rbuff2, this%gpC)
        call transpose_y_to_z(rbuff2,                   rbuff3, this%gpC)
        call this%compute_z_mean(rbuff3, this%u_mean)
        write(tempname,"(A3,I2.2,A5,I6.6,A6)") "Run",this%runID, "_um_t",this%step,".3Dstt"
        fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
        call decomp_2d_write_one(3, rbuff3, fname)

        ! v_avg
        call transpose_x_to_y(this%v_mean3D/tidSUMreal, rbuff2, this%gpC)
        call transpose_y_to_z(rbuff2,                   rbuff4, this%gpC)
        call this%compute_z_mean(rbuff4, this%v_mean)
        write(tempname,"(A3,I2.2,A5,I6.6,A6)") "Run",this%runID, "_vm_t",this%step,".3Dstt"
        fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
        call decomp_2d_write_one(3, rbuff4, fname)
      
        ! w_avg
        call transpose_x_to_y(this%w_mean3D/tidSUMreal, rbuff2, this%gpC)
        call transpose_y_to_z(rbuff2,                   rbuff5, this%gpC)
        call this%compute_z_mean(rbuff5, this%w_mean)
        write(tempname,"(A3,I2.2,A5,I6.6,A6)") "Run",this%runID, "_wm_t",this%step,".3Dstt"
        fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
        call decomp_2d_write_one(3, rbuff5, fname)
      
        ! uu_avg - u_avg*u_avg - Total (1D) and Reynolds (3D)
        call transpose_x_to_y(this%uu_mean3D/tidSUMreal, rbuff2, this%gpC)
        call transpose_y_to_z(rbuff2,                    rbuff6, this%gpC)
        call this%compute_z_mean(rbuff6, this%uu_mean)
        this%uu_mean = this%uu_mean - this%u_mean*this%u_mean    !--Total
        rbuff6 = rbuff6 - rbuff3 * rbuff3                        !--Reynolds
        write(tempname,"(A3,I2.2,A6,I6.6,A6)") "Run",this%runID, "_uum_t",this%step,".3Dstt"
        fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
        call decomp_2d_write_one(3, rbuff6, fname)
      
        ! uv_avg - u_avg*v_avg - Total (1D) and Reynolds (3D)
        call transpose_x_to_y(this%uv_mean3D/tidSUMreal, rbuff2, this%gpC)
        call transpose_y_to_z(rbuff2,                    rbuff6, this%gpC)
        call this%compute_z_mean(rbuff6, this%uv_mean)
        this%uv_mean = this%uv_mean - this%u_mean*this%v_mean    !--Total
        rbuff6 = rbuff6 - rbuff3 * rbuff4                        !--Reynolds
        write(tempname,"(A3,I2.2,A6,I6.6,A6)") "Run",this%runID, "_uvm_t",this%step,".3Dstt"
        fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
        call decomp_2d_write_one(3, rbuff6, fname)
      
        ! uw_avg - u_avg*w_avg - Reynolds
        call transpose_x_to_y(this%uw_mean3D/tidSUMreal, rbuff2, this%gpC)
        call transpose_y_to_z(rbuff2,                    rbuff6, this%gpC)
        rbuff6 = rbuff6 - rbuff3*rbuff5
        call this%compute_z_mean(rbuff6, this%uw_mean)
        write(tempname,"(A3,I2.2,A6,I6.6,A6)") "Run",this%runID, "_uwm_t",this%step,".3Dstt"
        fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
        call decomp_2d_write_one(3, rbuff6, fname)
      
        ! uw_avg - u_avg*w_avg - Dispersive
        rbuff6 = rbuff3*rbuff5
        call this%compute_z_mean(rbuff6, this%disperuw_mean)
        this%disperuw_mean = this%disperuw_mean - this%u_mean*this%w_mean
      
        ! vv_avg - v_avg*v_avg - Total (1D) and Reynolds (3D)
        call transpose_x_to_y(this%vv_mean3D/tidSUMreal, rbuff2, this%gpC)
        call transpose_y_to_z(rbuff2,                    rbuff6, this%gpC)
        call this%compute_z_mean(rbuff6, this%vv_mean)
        this%vv_mean = this%vv_mean - this%v_mean*this%v_mean    !--Total
        rbuff6 = rbuff6 - rbuff4 * rbuff4                        !--Reynolds
        write(tempname,"(A3,I2.2,A6,I6.6,A6)") "Run",this%runID, "_vvm_t",this%step,".3Dstt"
        fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
        call decomp_2d_write_one(3, rbuff6, fname)
      
        ! vw_avg - v_avg*w_avg - Reynolds
        call transpose_x_to_y(this%vw_mean3D/tidSUMreal, rbuff2, this%gpC)
        call transpose_y_to_z(rbuff2,                    rbuff6, this%gpC)
        rbuff6 = rbuff6 - rbuff4*rbuff5
        call this%compute_z_mean(rbuff6, this%vw_mean)
        write(tempname,"(A3,I2.2,A6,I6.6,A6)") "Run",this%runID, "_vwm_t",this%step,".3Dstt"
        fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
        call decomp_2d_write_one(3, rbuff6, fname)
      
        ! vw_avg - v_avg*w_avg - Dispersive
        rbuff6 = rbuff4*rbuff5
        call this%compute_z_mean(rbuff6, this%dispervw_mean)
        this%dispervw_mean = this%dispervw_mean - this%v_mean*this%w_mean
      
        ! ww_avg - w_avg*w_avg - Total (1D) and Reynolds (3D)
        call transpose_x_to_y(this%ww_mean3D/tidSUMreal, rbuff2, this%gpC)
        call transpose_y_to_z(rbuff2,                    rbuff6, this%gpC)
        call this%compute_z_mean(rbuff6, this%ww_mean)
        this%ww_mean = this%ww_mean - this%w_mean*this%w_mean    !--Total
        rbuff6 = rbuff6 - rbuff5 * rbuff5                        !--Reynolds
        write(tempname,"(A3,I2.2,A6,I6.6,A6)") "Run",this%runID, "_wwm_t",this%step,".3Dstt"
        fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
        call decomp_2d_write_one(3, rbuff6, fname)
      
        !-----Turbulent transport of TKE budget-------
          ! triple correlation for transport term in TKE budget
          ! x term in xdecomp
          !rbuff1 = this%u_mean3D*this%u_mean3D + this%v_mean3D*this%v_mean3D + this%w_mean3D*this%w_mean3D
          rbuff1 = this%tketurbtranspx_mean3D/tidSUMreal - two*(this%u_mean3D*this%uu_mean3D + this%v_mean3D*this%uv_mean3D + this%w_mean3D*this%uw_mean3D)/(tidSumreal**2) &
                 + ( two*(this%u_mean3D*this%u_mean3D + this%v_mean3D*this%v_mean3D + this%w_mean3D*this%w_mean3D)/tidSUMreal**2 - (this%uu_mean3D + this%vv_mean3D + this%ww_mean3D)/tidSUMreal )*(this%u_mean3D/tidSUMreal)
          call this%spectC%fft(rbuff1,cbuffy1)   
          call this%spectC%mTimes_ik1_ip(cbuffy1)
          call this%spectC%ifft(cbuffy1,rbuff0)

          ! add y term in x decomp
          rbuff1 = this%tketurbtranspy_mean3D/tidSUMreal - two*(this%u_mean3D*this%uv_mean3D + this%v_mean3D*this%vv_mean3D + this%w_mean3D*this%vw_mean3D)/(tidSumreal**2) &
                 + ( two*(this%u_mean3D*this%u_mean3D + this%v_mean3D*this%v_mean3D + this%w_mean3D*this%w_mean3D)/tidSUMreal**2 - (this%uu_mean3D + this%vv_mean3D + this%ww_mean3D)/tidSUMreal )*(this%v_mean3D/tidSUMreal)
          call this%spectC%fft(rbuff1,cbuffy1)   
          call this%spectC%mTimes_ik2_ip(cbuffy1)
          call this%spectC%ifft(cbuffy1,rbuff1)
          rbuff0 = rbuff0 + rbuff1

          ! take sum of x and y terms to z decomp
          call transpose_x_to_y(rbuff0, rbuff2, this%gpC)
          call transpose_y_to_z(rbuff2, rbuff4, this%gpC)

          ! compute z term in z decomp
          rbuff1 = this%tketurbtranspz_mean3D/tidSUMreal - two*(this%u_mean3D*this%uw_mean3D + this%v_mean3D*this%vw_mean3D + this%w_mean3D*this%ww_mean3D)/(tidSumreal**2) &
                 + ( two*(this%u_mean3D*this%u_mean3D + this%v_mean3D*this%v_mean3D + this%w_mean3D*this%w_mean3D)/tidSUMreal**2 - (this%uu_mean3D + this%vv_mean3D + this%ww_mean3D)/tidSUMreal )*(this%w_mean3D/tidSUMreal)
          call transpose_x_to_y(rbuff1, rbuff2, this%gpC)
          call transpose_y_to_z(rbuff2, rbuff3, this%gpC)
          ! interpolate rbuff3 from C to E
          call this%Pade6opZ%interpz_C2E(rbuff3, rbuff3E, 0,0)
          call this%Pade6opZ%ddz_E2C(rbuff3E,rbuff3,0,0)

          ! add x and y terms to z term
          rbuff3 = rbuff3 + rbuff4
          rbuff3 = -half*rbuff3

          call this%compute_z_mean(rbuff3, this%tkett_mean)
          write(tempname,"(A3,I2.2,A7,I6.6,A6)") "Run",this%runID, "_tktt_t",this%step,".3Dstt"
          fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
          call decomp_2d_write_one(3, rbuff3, fname)
        !-----Done turbulent transport of TKE budget-------

        !-----Mean advection term of TKE budget------
          ! compute tke first
          rbuff1 = this%u_mean3D*this%u_mean3D + this%v_mean3D*this%v_mean3D + this%w_mean3D*this%w_mean3D
          rbuff1 = -rbuff1/tidSUMreal**2
          rbuff1 = rbuff1 + (this%uu_mean3D + this%vv_mean3D + this%ww_mean3D)/tidSUMreal

          ! transpose to z and take z derivative
          call transpose_x_to_y(rbuff1, rbuff2, this%gpC)
          call transpose_y_to_z(rbuff2, rbuff3, this%gpC)
          call this%Pade6opZ%ddz_C2E(rbuff3,rbuff3E,0,0)
          
          ! transpose w_mean3D to z, interpolate to E and multiply
          call transpose_x_to_y(this%w_mean3D/tidSUMreal, rbuff2, this%gpC)
          call transpose_y_to_z(rbuff2,                   rbuff4, this%gpC)
          call this%Pade6opZ%interpz_C2E(rbuff4,rbuff4E,0,0)
          rbuff4E = rbuff4E * rbuff3E

          ! interpolate E to C
          call this%Pade6opZ%interpz_E2C(rbuff4E, rbuff4, 0,0)


          ! x derivative
          call this%spectC%fft(rbuff1,cbuffy1)   
          call this%spectC%mTimes_ik1_oop(cbuffy1, cbuffy2)
          call this%spectC%ifft(cbuffy2,rbuff0)

          ! y derivative
          call this%spectC%fft(rbuff1,cbuffy1)   
          call this%spectC%mTimes_ik2_oop(cbuffy1, cbuffy2)
          call this%spectC%ifft(cbuffy2,rbuff1)
          rbuff0 = this%u_mean3D*rbuff0 + this%v_mean3D*rbuff1

          ! transpose sum of x and y parts to z and add z part
          call transpose_x_to_y(rbuff0/tidSUMreal, rbuff2, this%gpC)
          call transpose_y_to_z(rbuff2,            rbuff3, this%gpC)
          rbuff3 = rbuff3 + rbuff4
          rbuff3 = -half*rbuff3

          ! write outputs
          call this%compute_z_mean(rbuff3, this%tkeadv_mean)
          write(tempname,"(A3,I2.2,A7,I6.6,A6)") "Run",this%runID, "_tkma_t",this%step,".3Dstt"
          fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
          call decomp_2d_write_one(3, rbuff3, fname)
        !-----Done mean advection term of TKE budget------

        S11_mean3D => Stmp(:,:,:,1);   S12_mean3D => Stmp(:,:,:,2);    S13_mean3D => Stmp(:,:,:,3) 
                                       S22_mean3D => Stmp(:,:,:,4);    S23_mean3D => Stmp(:,:,:,5) 
                                                                       S33_mean3D => Stmp(:,:,:,6) 
        call this%compute_Sijmean(Stmp)

        ! ---- Shear production in TKE budget equation-------
          rbuff1 = - (this%uu_mean3D/tidSUMreal - this%u_mean3D*this%u_mean3D/tidSUMreal**2)*S11_mean3D &
                   - (this%vv_mean3D/tidSUMreal - this%v_mean3D*this%v_mean3D/tidSUMreal**2)*S22_mean3D &
                   - (this%ww_mean3D/tidSUMreal - this%w_mean3D*this%w_mean3D/tidSUMreal**2)*S33_mean3D
          rbuff1 = rbuff1 - two*( &
                   + (this%uv_mean3D/tidSUMreal - this%u_mean3D*this%v_mean3D/tidSUMreal**2)*S12_mean3D &
                   + (this%uw_mean3D/tidSUMreal - this%u_mean3D*this%w_mean3D/tidSUMreal**2)*S13_mean3D &
                   + (this%vw_mean3D/tidSUMreal - this%v_mean3D*this%w_mean3D/tidSUMreal**2)*S23_mean3D )

          call transpose_x_to_y(rbuff1, rbuff2, this%gpC)
          call transpose_y_to_z(rbuff2, rbuff3, this%gpC)
          call this%compute_z_mean(rbuff3, this%tkeprod_mean)
          write(tempname,"(A3,I2.2,A7,I6.6,A6)") "Run",this%runID, "_tkpr_t",this%step,".3Dstt"
          fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
          call decomp_2d_write_one(3, rbuff3, fname)
        ! ---- Done Shear production in TKE budget equation-------
      
        ! ---- Turbulent dissipation in MKE budget equation-------
          ! --- Identical to negative of shear production in TKE equation---
          this%mkedisp_mean = -this%tkeprod_mean
        ! ---- Done Turbulent dissipation in MKE budget equation-------
        
        !-----Transport of turbulent stresses in MKE budget-------
          ! x term in xdecomp
          rbuff1 = - (this%uu_mean3D/tidSUMreal - this%u_mean3D*this%u_mean3D/tidSUMreal**2)*this%u_mean3D/tidSUMreal &
                   - (this%uv_mean3D/tidSUMreal - this%u_mean3D*this%v_mean3D/tidSUMreal**2)*this%v_mean3D/tidSUMreal &
                   - (this%uw_mean3D/tidSUMreal - this%u_mean3D*this%w_mean3D/tidSUMreal**2)*this%w_mean3D/tidSUMreal
          call this%spectC%fft(rbuff1,cbuffy1)   
          call this%spectC%mTimes_ik1_ip(cbuffy1)
          call this%spectC%ifft(cbuffy1,rbuff0)

          ! y term in xdecomp
          rbuff1 = - (this%uv_mean3D/tidSUMreal - this%u_mean3D*this%v_mean3D/tidSUMreal**2)*this%u_mean3D/tidSUMreal &
                   - (this%vv_mean3D/tidSUMreal - this%v_mean3D*this%v_mean3D/tidSUMreal**2)*this%v_mean3D/tidSUMreal &
                   - (this%vw_mean3D/tidSUMreal - this%v_mean3D*this%w_mean3D/tidSUMreal**2)*this%w_mean3D/tidSUMreal
          call this%spectC%fft(rbuff1,cbuffy1)   
          call this%spectC%mTimes_ik2_oop(cbuffy1, cbuffy2)
          call this%spectC%ifft(cbuffy2,rbuff1)

          ! transpose sum of x and y parts to z
          call transpose_x_to_y(rbuff0+rbuff1, rbuff2, this%gpC)
          call transpose_y_to_z(rbuff2,        rbuff4, this%gpC)

          ! compute z term in z decomp
          rbuff1 = - (this%uw_mean3D/tidSUMreal - this%u_mean3D*this%w_mean3D/tidSUMreal**2)*this%u_mean3D/tidSUMreal &
                   - (this%vw_mean3D/tidSUMreal - this%v_mean3D*this%w_mean3D/tidSUMreal**2)*this%v_mean3D/tidSUMreal &
                   - (this%ww_mean3D/tidSUMreal - this%w_mean3D*this%w_mean3D/tidSUMreal**2)*this%w_mean3D/tidSUMreal
          call transpose_x_to_y(rbuff1, rbuff2, this%gpC)
          call transpose_y_to_z(rbuff2, rbuff3, this%gpC)
          ! interpolate rbuff3 from C to E
          call this%Pade6opZ%interpz_C2E(rbuff3, rbuff3E, 0,0)
          call this%Pade6opZ%ddz_E2C(rbuff3E,rbuff3,0,0)

          ! add x and y terms to z term
          rbuff3 = rbuff3 + rbuff4

          call this%compute_z_mean(rbuff3, this%mkett_mean)
          write(tempname,"(A3,I2.2,A7,I6.6,A6)") "Run",this%runID, "_mktt_t",this%step,".3Dstt"
          fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
          call decomp_2d_write_one(3, rbuff3, fname)
        !-----Done Transport of turbulent stresses in MKE budget-------

        !-----Mean advection term of MKE budget------
          ! compute mke first
          rbuff1 = this%u_mean3D*this%u_mean3D + this%v_mean3D*this%v_mean3D + this%w_mean3D*this%w_mean3D
          rbuff1 = rbuff1/tidSUMreal**2

          ! transpose to z and take z derivative
          call transpose_x_to_y(rbuff1, rbuff2, this%gpC)
          call transpose_y_to_z(rbuff2, rbuff3, this%gpC)
          call this%Pade6opZ%ddz_C2E(rbuff3,rbuff3E,0,0)
          
          ! transpose w_mean3D to z, interpolate to E and multiply
          call transpose_x_to_y(this%w_mean3D/tidSUMreal, rbuff2, this%gpC)
          call transpose_y_to_z(rbuff2,                   rbuff4, this%gpC)
          call this%Pade6opZ%interpz_C2E(rbuff4,rbuff4E,0,0)
          rbuff4E = rbuff4E * rbuff3E

          ! interpolate E to C
          call this%Pade6opZ%interpz_E2C(rbuff4E, rbuff4, 0,0)


          ! x derivative
          call this%spectC%fft(rbuff1,cbuffy1)   
          call this%spectC%mTimes_ik1_oop(cbuffy1, cbuffy2)
          call this%spectC%ifft(cbuffy2,rbuff0)

          ! y derivative
          call this%spectC%fft(rbuff1,cbuffy1)   
          call this%spectC%mTimes_ik2_oop(cbuffy1, cbuffy2)
          call this%spectC%ifft(cbuffy2,rbuff1)
          rbuff0 = this%u_mean3D*rbuff0 + this%v_mean3D*rbuff1

          ! transpose sum of x and y parts to z and add z part
          call transpose_x_to_y(rbuff0/tidSUMreal, rbuff2, this%gpC)
          call transpose_y_to_z(rbuff2,            rbuff3, this%gpC)
          rbuff3 = rbuff3 + rbuff4
          rbuff3 = -half*rbuff3

          ! write outputs
          call this%compute_z_mean(rbuff3, this%mkeadv_mean)
          write(tempname,"(A3,I2.2,A7,I6.6,A6)") "Run",this%runID, "_mkma_t",this%step,".3Dstt"
          fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
          call decomp_2d_write_one(3, rbuff3, fname)
        !-----Done mean advection term of MKE budget------

        if(this%fastCalcPressure .or. this%storePressure) then
            ! p_avg
            call transpose_x_to_y(this%p_mean3D/tidSUMreal, rbuff2, this%gpC)
            call transpose_y_to_z(rbuff2,                   rbuff6, this%gpC)
            call this%compute_z_mean(rbuff6, this%p_mean)
            write(tempname,"(A3,I2.2,A5,I6.6,A6)") "Run",this%runID, "_pm_t",this%step,".3Dstt"
            fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
            call decomp_2d_write_one(3, rbuff6, fname)

            ! compute ddxj(p_avg * uavg_j) -- Pressure transport term in mean KE eqn
              ! ddx(p_avg*uavg_1) in x decomp
              rbuff1 = this%p_mean3D*this%u_mean3D/(tidSUMreal**2)
              call this%spectC%fft(rbuff1,cbuffy1)   
              call this%spectC%mTimes_ik1_oop(cbuffy1, cbuffy2)
              call this%spectC%ifft(cbuffy2,rbuff0)

              ! ddy(p_avg*uavg_2) in x decomp
              rbuff1 = this%p_mean3D*this%v_mean3D/(tidSUMreal**2)
              call this%spectC%fft(rbuff1,cbuffy1)   
              call this%spectC%mTimes_ik2_oop(cbuffy1, cbuffy2)
              call this%spectC%ifft(cbuffy2,rbuff1)
              rbuff0 = rbuff0 + rbuff1

              ! take sum of x and y terms to z decomp
              call transpose_x_to_y(rbuff0, rbuff2, this%gpC)
              call transpose_y_to_z(rbuff2, rbuff4, this%gpC)

              ! ddz(p_avg*uavg_3) in z decomp
              rbuff1 = this%p_mean3D*this%w_mean3D/(tidSUMreal**2)
              ! transpose to z and take z derivative
              call transpose_x_to_y(rbuff1, rbuff2, this%gpC)
              call transpose_y_to_z(rbuff2, rbuff3, this%gpC)
              call this%Pade6opZ%ddz_C2E(rbuff3,rbuff3E,0,0)
              call this%Pade6opZ%interpz_E2C(rbuff3E, rbuff3, 0,0)       

              rbuff3 = -(rbuff3 + rbuff4)

            call this%compute_z_mean(rbuff3, this%mkept_mean)       ! Pressure transport mean (term in mean KE eqn)
            write(tempname,"(A3,I2.2,A7,I6.6,A6)") "Run",this%runID, "_mkpt_t",this%step,".3Dstt"
            fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
            call decomp_2d_write_one(3, rbuff3, fname)
            ! Done computing ddxj(p_avg * uavg_j) -- Pressure transport term in mean KE eqn

      
            ! compute ddxj(p' * u'_j) -- Pressure transport term in turbulent KE eqn
              ! pu_avg - u_avg*p_avg
              rbuff1 = this%pu_mean3D/tidSUMreal - this%u_mean3D * this%p_mean3D / tidSUMreal**2
              call this%spectC%fft(rbuff1,cbuffy1)   
              call this%spectC%mTimes_ik1_oop(cbuffy1, cbuffy2)
              call this%spectC%ifft(cbuffy2,rbuff0)

              ! pv_avg - v_avg*p_avg - Reynolds only
              rbuff1 = this%pv_mean3D/tidSUMreal - this%v_mean3D * this%p_mean3D / tidSUMreal**2
              call this%spectC%fft(rbuff1,cbuffy1)   
              call this%spectC%mTimes_ik2_oop(cbuffy1, cbuffy2)
              call this%spectC%ifft(cbuffy2,rbuff1)
              rbuff0 = rbuff0 + rbuff1

              ! take sum of x and y terms to z decomp
              call transpose_x_to_y(rbuff0, rbuff2, this%gpC)
              call transpose_y_to_z(rbuff2, rbuff4, this%gpC)

              ! pw_avg - w_avg*p_avg - Reynolds only
              rbuff1 = this%pw_mean3D/tidSUMreal - this%w_mean3D * this%p_mean3D / tidSUMreal**2
              call transpose_x_to_y(rbuff1, rbuff2, this%gpC)
              call transpose_y_to_z(rbuff2, rbuff3, this%gpC)
              call this%Pade6opZ%ddz_C2E(rbuff3,rbuff3E,0,0)
              call this%Pade6opZ%interpz_E2C(rbuff3E, rbuff3, 0,0)       

              rbuff3 = -(rbuff3 + rbuff4)

            call this%compute_z_mean(rbuff3, this%tkept_mean)       ! Pressure transport turbulent (term in turbulent KE eqn)
            write(tempname,"(A3,I2.2,A7,I6.6,A6)") "Run",this%runID, "_tkpt_t",this%step,".3Dstt"
            fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
            call decomp_2d_write_one(3, rbuff3, fname)
            ! Done computing ddxj(p' * u'_j) -- Pressure transport term in turbulent KE eqn

        endif

        if(.not. this%isInviscid) then
           ! mkevdif_mean; mkevdsp_mean; tkevdif_mean; tkevdsp_mean to be written
        endif

        if(this%useSGS) then
            ! tau11SGS_avg
            call transpose_x_to_y(this%tau11_mean3D/tidSUMreal, rbuff2, this%gpC)
            call transpose_y_to_z(rbuff2,                       rbuff3, this%gpC)
            call this%compute_z_mean(rbuff3, this%tau11_mean)
            write(tempname,"(A3,I2.2,A6,I6.6,A6)") "Run",this%runID, "_t11_t",this%step,".3Dstt"
            fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
            call decomp_2d_write_one(3, rbuff3, fname)
      
            ! tau12SGS_avg
            call transpose_x_to_y(this%tau12_mean3D/tidSUMreal, rbuff2, this%gpC)
            call transpose_y_to_z(rbuff2,                       rbuff3, this%gpC)
            call this%compute_z_mean(rbuff3, this%tau12_mean)
            write(tempname,"(A3,I2.2,A6,I6.6,A6)") "Run",this%runID, "_t12_t",this%step,".3Dstt"
            fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
            call decomp_2d_write_one(3, rbuff3, fname)
      
            ! tau13SGS_avg
            call transpose_x_to_y(this%tau13_mean3D/tidSUMreal, rbuff2, this%gpC)
            call transpose_y_to_z(rbuff2,                       rbuff3, this%gpC)
            call this%compute_z_mean(rbuff3, this%tau13_mean)
            write(tempname,"(A3,I2.2,A6,I6.6,A6)") "Run",this%runID, "_t13_t",this%step,".3Dstt"
            fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
            call decomp_2d_write_one(3, rbuff3, fname)
      
            ! tau22SGS_avg
            call transpose_x_to_y(this%tau22_mean3D/tidSUMreal, rbuff2, this%gpC)
            call transpose_y_to_z(rbuff2,                       rbuff3, this%gpC)
            call this%compute_z_mean(rbuff3, this%tau22_mean)
            write(tempname,"(A3,I2.2,A6,I6.6,A6)") "Run",this%runID, "_t22_t",this%step,".3Dstt"
            fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
            call decomp_2d_write_one(3, rbuff3, fname)
      
            ! tau23SGS_avg
            call transpose_x_to_y(this%tau23_mean3D/tidSUMreal, rbuff2, this%gpC)
            call transpose_y_to_z(rbuff2,                       rbuff3, this%gpC)
            call this%compute_z_mean(rbuff3, this%tau23_mean)
            write(tempname,"(A3,I2.2,A6,I6.6,A6)") "Run",this%runID, "_t23_t",this%step,".3Dstt"
            fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
            call decomp_2d_write_one(3, rbuff3, fname)
      
            ! tau33SGS_avg
            call transpose_x_to_y(this%tau33_mean3D/tidSUMreal, rbuff2, this%gpC)
            call transpose_y_to_z(rbuff2,                       rbuff3, this%gpC)
            call this%compute_z_mean(rbuff3, this%tau33_mean)
            write(tempname,"(A3,I2.2,A6,I6.6,A6)") "Run",this%runID, "_t33_t",this%step,".3Dstt"
            fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
            call decomp_2d_write_one(3, rbuff3, fname)
      
            ! SGS dissipation in MKE equation
            rbuff1 = (    this%tau11_mean3D*S11_mean3D + this%tau22_mean3D*S22_mean3D + this%tau33_mean3D*S33_mean3D + & 
                     two*(this%tau12_mean3D*S12_mean3D + this%tau13_mean3D*S13_mean3D + this%tau23_mean3D*S23_mean3D)  )/tidSUMreal
            call transpose_x_to_y(rbuff1, rbuff2, this%gpC)
            call transpose_y_to_z(rbuff2, rbuff3, this%gpC)
            call this%compute_z_mean(rbuff3, this%mkesgsd_mean)
            write(tempname,"(A3,I2.2,A9,I6.6,A6)") "Run",this%runID, "_mksgsd_t",this%step,".3Dstt"
            fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
            call decomp_2d_write_one(3, rbuff3, fname)

            ! SGS dissipation in TKE equation
            rbuff1 = this%sgsdissp_mean3D/tidSUMreal - rbuff1
            call transpose_x_to_y(rbuff1, rbuff2, this%gpC)
            call transpose_y_to_z(rbuff2, rbuff3, this%gpC)
            call this%compute_z_mean(rbuff3, this%tkesgsd_mean)
            write(tempname,"(A3,I2.2,A9,I6.6,A6)") "Run",this%runID, "_tksgsd_t",this%step,".3Dstt"
            fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
            call decomp_2d_write_one(3, rbuff3, fname)
      
            ! Mean transport of SGS stress in MKE equation
              ! x term in xdecomp
              rbuff1 = (this%u_mean3D*this%tau11_mean3D + this%v_mean3D*this%tau12_mean3D + this%w_mean3D*this%tau13_mean3D)/tidSUMreal**2
              call this%spectC%fft(rbuff1,cbuffy1)   
              call this%spectC%mTimes_ik1_oop(cbuffy1, cbuffy2)
              call this%spectC%ifft(cbuffy2,rbuff0)

              ! y term in xdecomp
              rbuff1 = (this%u_mean3D*this%tau12_mean3D + this%v_mean3D*this%tau22_mean3D + this%w_mean3D*this%tau23_mean3D)/tidSUMreal**2
              call this%spectC%fft(rbuff1,cbuffy1)   
              call this%spectC%mTimes_ik2_oop(cbuffy1, cbuffy2)
              call this%spectC%ifft(cbuffy2,rbuff1)
              rbuff0 = rbuff0 + rbuff1

              ! take sum of x and y terms to z decomp
              call transpose_x_to_y(rbuff0, rbuff2, this%gpC)
              call transpose_y_to_z(rbuff2, rbuff4, this%gpC)

              ! z term in zdecomp
              rbuff1 = (this%u_mean3D*this%tau13_mean3D + this%v_mean3D*this%tau23_mean3D + this%w_mean3D*this%tau33_mean3D)/tidSUMreal**2
              call transpose_x_to_y(rbuff1, rbuff2, this%gpC)
              call transpose_y_to_z(rbuff2, rbuff3, this%gpC)
              call this%Pade6opZ%ddz_C2E(rbuff3,rbuff3E,0,0)
              call this%Pade6opZ%interpz_E2C(rbuff3E, rbuff3, 0,0)       

              rbuff3 = -(rbuff3 + rbuff4)

              call this%compute_z_mean(rbuff3, this%mkesgst_mean)       ! Pressure transport turbulent (term in turbulent KE eqn)
              write(tempname,"(A3,I2.2,A9,I6.6,A6)") "Run",this%runID, "_mksgst_t",this%step,".3Dstt"
              fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
              call decomp_2d_write_one(3, rbuff3, fname)
            ! Done mean transport of SGS stress in MKE equation

            ! Turbulent transport of SGS stress in TKE equation
              ! x term in xdecomp
              rbuff1 = this%tauu1_mean3D/tidSUMreal - (this%u_mean3D*this%tau11_mean3D + this%v_mean3D*this%tau12_mean3D + this%w_mean3D*this%tau13_mean3D)/tidSUMreal**2
              call this%spectC%fft(rbuff1,cbuffy1)   
              call this%spectC%mTimes_ik1_oop(cbuffy1, cbuffy2)
              call this%spectC%ifft(cbuffy2,rbuff0)

              ! y term in xdecomp
              rbuff1 = this%tauu2_mean3D/tidSUMreal - (this%u_mean3D*this%tau12_mean3D + this%v_mean3D*this%tau22_mean3D + this%w_mean3D*this%tau23_mean3D)/tidSUMreal**2
              call this%spectC%fft(rbuff1,cbuffy1)   
              call this%spectC%mTimes_ik2_oop(cbuffy1, cbuffy2)
              call this%spectC%ifft(cbuffy2,rbuff1)
              rbuff0 = rbuff0 + rbuff1

              ! take sum of x and y terms to z decomp
              call transpose_x_to_y(rbuff0, rbuff2, this%gpC)
              call transpose_y_to_z(rbuff2, rbuff4, this%gpC)

              ! z term in zdecomp
              rbuff1 = this%tauu3_mean3D/tidSUMreal - (this%u_mean3D*this%tau13_mean3D + this%v_mean3D*this%tau23_mean3D + this%w_mean3D*this%tau33_mean3D)/tidSUMreal**2
              call transpose_x_to_y(rbuff1, rbuff2, this%gpC)
              call transpose_y_to_z(rbuff2, rbuff3, this%gpC)
              call this%Pade6opZ%ddz_C2E(rbuff3,rbuff3E,0,0)
              call this%Pade6opZ%interpz_E2C(rbuff3E, rbuff3, 0,0)       

              rbuff3 = -(rbuff3 + rbuff4)

              call this%compute_z_mean(rbuff3, this%tkesgst_mean)       ! Pressure transport turbulent (term in turbulent KE eqn)
              write(tempname,"(A3,I2.2,A9,I6.6,A6)") "Run",this%runID, "_tksgst_t",this%step,".3Dstt"
              fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
              call decomp_2d_write_one(3, rbuff3, fname)
            ! Done turbulent transport of SGS stress in TKE equation
        endif

        if(this%useWindTurbines) then
           !----Turbine work term in MKE budget-------
           rbuff1 = ( this%u_mean3D*this%turbfx_mean3D + this%v_mean3D*this%turbfy_mean3D + this%w_mean3D*this%turbfy_mean3D ) / tidSUMreal**2
           call transpose_x_to_y(rbuff1, rbuff2, this%gpC)
           call transpose_y_to_z(rbuff2, rbuff3, this%gpC)
           call this%compute_z_mean(rbuff3, this%mketurbf_mean)
           write(tempname,"(A3,I2.2,A9,I6.6,A6)") "Run",this%runID, "_mktrbf_t",this%step,".3Dstt"
           fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
           call decomp_2d_write_one(3, rbuff3, fname)
           !----Done turbine work term in MKE budget-------

           !----Turbine work term in TKE budget-------
           rbuff1 = this%uturbf_mean3D/tidSUMreal -  ( this%u_mean3D*this%turbfx_mean3D + &
                    this%v_mean3D*this%turbfy_mean3D + this%w_mean3D*this%turbfy_mean3D ) / tidSUMreal**2
           call transpose_x_to_y(rbuff1, rbuff2, this%gpC)
           call transpose_y_to_z(rbuff2, rbuff3, this%gpC)
           call this%compute_z_mean(rbuff3, this%tketurbf_mean)
           write(tempname,"(A3,I2.2,A9,I6.6,A6)") "Run",this%runID, "_tktrbf_t",this%step,".3Dstt"
           fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
           call decomp_2d_write_one(3, rbuff3, fname)
           !----Done turbine work term in TKE budget-------

           ! write out turbfx_mean3D
           call transpose_x_to_y(this%turbfx_mean3D/tidSUMreal, rbuff2, this%gpC)
           call transpose_y_to_z(rbuff2,                        rbuff3, this%gpC)
           call this%compute_z_mean(rbuff3, this%turbfx_mean)
           write(tempname,"(A3,I2.2,A7,I6.6,A6)") "Run",this%runID, "_trbx_t",this%step,".3Dstt"
           fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
           call decomp_2d_write_one(3, rbuff3, fname)

           ! write out turbfy_mean3D
           call transpose_x_to_y(this%turbfy_mean3D/tidSUMreal, rbuff2, this%gpC)
           call transpose_y_to_z(rbuff2,                        rbuff3, this%gpC)
           call this%compute_z_mean(rbuff3, this%turbfy_mean)
           write(tempname,"(A3,I2.2,A7,I6.6,A6)") "Run",this%runID, "_trby_t",this%step,".3Dstt"
           fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
           call decomp_2d_write_one(3, rbuff3, fname)

           ! write out turbfz_mean3D
           call transpose_x_to_y(this%turbfz_mean3D/tidSUMreal, rbuff2, this%gpC)
           call transpose_y_to_z(rbuff2,                        rbuff3, this%gpC)
           call this%compute_z_mean(rbuff3, this%turbfz_mean)
           write(tempname,"(A3,I2.2,A7,I6.6,A6)") "Run",this%runID, "_trbz_t",this%step,".3Dstt"
           fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
           call decomp_2d_write_one(3, rbuff3, fname)
        endif

        if(this%isStratified) then
            ! T_avg
            call transpose_x_to_y(this%T_mean3D/tidSUMreal, rbuff2, this%gpC)
            call transpose_y_to_z(rbuff2,                   rbuff6, this%gpC)
            call this%compute_z_mean(rbuff6, this%T_mean)
            write(tempname,"(A3,I2.2,A5,I6.6,A6)") "Run",this%runID, "_Tm_t",this%step,".3Dstt"
            fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
            call decomp_2d_write_one(3, rbuff6, fname)
      
            ! uT_avg - u_avg*T_avg - Reynolds only
            rbuff1 = this%uT_mean3D/tidSUMreal - this%u_mean3D * this%T_mean3D / tidSUMreal**2
            call transpose_x_to_y(rbuff1, rbuff2, this%gpC)
            call transpose_y_to_z(rbuff2, rbuff3, this%gpC)
            call this%compute_z_mean(rbuff3, this%uT_mean)       !--Reynolds
            write(tempname,"(A3,I2.2,A6,I6.6,A6)") "Run",this%runID, "_uTm_t",this%step,".3Dstt"
            fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
            call decomp_2d_write_one(3, rbuff3, fname)
      
            ! vT_avg - v_avg*T_avg - Reynolds only
            rbuff1 = this%vT_mean3D/tidSUMreal - this%v_mean3D * this%T_mean3D / tidSUMreal**2
            call transpose_x_to_y(rbuff1, rbuff2, this%gpC)
            call transpose_y_to_z(rbuff2, rbuff3, this%gpC)
            call this%compute_z_mean(rbuff3, this%vT_mean)
            write(tempname,"(A3,I2.2,A6,I6.6,A6)") "Run",this%runID, "_vTm_t",this%step,".3Dstt"
            fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
            call decomp_2d_write_one(3, rbuff3, fname)
      
            ! wT_avg - w_avg*T_avg - Reynolds only
            rbuff1 = this%wT_mean3D/tidSUMreal - this%w_mean3D * this%T_mean3D / tidSUMreal**2
            call transpose_x_to_y(rbuff1, rbuff2, this%gpC)
            call transpose_y_to_z(rbuff2, rbuff3, this%gpC)
            call this%compute_z_mean(rbuff3, this%wT_mean)
            write(tempname,"(A3,I2.2,A6,I6.6,A6)") "Run",this%runID, "_wTm_t",this%step,".3Dstt"
            fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
            call decomp_2d_write_one(3, rbuff3, fname)
      
            ! TT_avg - T_avg*T_avg - Reynolds only
            rbuff1 = this%TT_mean3D/tidSUMreal - this%T_mean3D * this%T_mean3D / tidSUMreal**2
            call transpose_x_to_y(rbuff1, rbuff2, this%gpC)
            call transpose_y_to_z(rbuff2, rbuff3, this%gpC)
            call this%compute_z_mean(rbuff3, this%TT_mean)
            write(tempname,"(A3,I2.2,A6,I6.6,A6)") "Run",this%runID, "_TTm_t",this%step,".3Dstt"
            fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
            call decomp_2d_write_one(3, rbuff3, fname)

            if(this%isStratified) then
                ! q1_mean
                call transpose_x_to_y(this%q1_mean3D/tidSUMreal, rbuff2, this%gpC)
                call transpose_y_to_z(rbuff2,                    rbuff3, this%gpC)
                call this%compute_z_mean(rbuff3, this%q1_mean)
                write(tempname,"(A3,I2.2,A5,I6.6,A6)") "Run",this%runID, "_q1_t",this%step,".3Dstt"
                fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                call decomp_2d_write_one(3, rbuff3, fname)
     
                ! q2_mean
                call transpose_x_to_y(this%q2_mean3D/tidSUMreal, rbuff2, this%gpC)
                call transpose_y_to_z(rbuff2,                    rbuff3, this%gpC)
                call this%compute_z_mean(rbuff3, this%q2_mean)
                write(tempname,"(A3,I2.2,A5,I6.6,A6)") "Run",this%runID, "_q2_t",this%step,".3Dstt"
                fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                call decomp_2d_write_one(3, rbuff3, fname)
     
                ! q3_mean
                call transpose_x_to_y(this%q3_mean3D/tidSUMreal, rbuff2, this%gpC)
                call transpose_y_to_z(rbuff2,                    rbuff3, this%gpC)
                call this%compute_z_mean(rbuff3, this%q3_mean)
                write(tempname,"(A3,I2.2,A5,I6.6,A6)") "Run",this%runID, "_q3_t",this%step,".3Dstt"
                fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                call decomp_2d_write_one(3, rbuff3, fname)
            endif
        endif

        call message(1, "Dumped 3D stats files")

        nullify(rbuff1,rbuff2,rbuff3,rbuff4,rbuff5,rbuff6, rbuff0, cbuffy1, cbuffy2, rbuff3E, rbuff4E)

        ! dump horizontal averages
        if(this%useWindTurbines) then
            this%runningSum_turb = zero
            call MPI_reduce(this%runningSum_sc_turb, this%runningSum_turb, 8*this%WindTurbineArr%nTurbines, mpirkind, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        endif
        if (nrank == 0) then
            write(tempname,"(A3,I2.2,A2,I6.6,A4)") "Run", this%RunID,"_t",tid,".stt"
            fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
            call write_2d_ascii(this%horzavgstats,fname)

            write(tempname,"(A3,I2.2,A2,I6.6,A4)") "Run", this%RunID,"_t",tid,".sth"   ! time and horz averages of scalars
            fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
            open(unit=771,file=fname,status='unknown')
            if(this%useWindTurbines) then
                write(771,'(e19.12,1x,i7,1x,8008(e19.12,1x))') this%tsim, this%tidSUM, this%runningSum_sc/tidSUMreal, this%runningSum_turb/tidSUMreal ! change if using more than 1000 turbines
            else
                write(771,'(e19.12,1x,i7,1x,5(e19.12,1x))') this%tsim, this%tidSUM, this%runningSum_sc/tidSUMreal
            endif
            close(771)
        end if
        !write(200+nrank,*) 4, this%tsim, this%runningSum_sc(2)/tidSUMreal
        call message(1, "Just dumped a .stt file")
        call message(2, "Number ot tsteps averaged:",this%tidSUM)

        if (this%computeSpectra) then
            ! Dump horizontally averaged x-spectra
            dirid = 2; decompdir = 2

            ! --- only 4, 5 or 6 planes of xspextra_mean in y-direction are being used
            nspectra = 4
            if(this%isStratified)                             nspectra = nspectra + 1
            if(this%fastCalcPressure .or. this%storePressure) nspectra = nspectra+1

            ! --- for k1 = 1, multiplication factor is 1.0,
            ! --- for k1 = 2:Nx/2+1, multiplication factor is 2.0
            normfac = two/real(size(this%cbuffyC(:,:,:,1),2),rkind)/tidSUMreal
            tmpvar(1:this%sp_gpC%ysz(1),1:nspectra,:) = normfac*this%xspectra_mean(1:this%sp_gpC%ysz(1),1:nspectra,:)
            if(this%sp_gpC%yst(1)==1) then
                tmpvar(1,1:nspectra,:) = half*tmpvar(1,1:nspectra,:)
            endif

            jindx = 1 ! u
            write(tempname,"(A3,I2.2,A8,I6.6,A4)") "Run", this%RunID,"_specu_t",tid,".out"
            fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
            call decomp_2d_write_plane(decompdir, tmpvar, dirid, jindx, fname, this%sp_gpC)

            jindx = 2 ! v
            write(tempname,"(A3,I2.2,A8,I6.6,A4)") "Run", this%RunID,"_specv_t",tid,".out"
            fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
            call decomp_2d_write_plane(decompdir, tmpvar, dirid, jindx, fname, this%sp_gpC)

            jindx = 3 ! w
            write(tempname,"(A3,I2.2,A8,I6.6,A4)") "Run", this%RunID,"_specw_t",tid,".out"
            fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
            call decomp_2d_write_plane(decompdir, tmpvar, dirid, jindx, fname, this%sp_gpC)

            jindx = 4 ! KE
            write(tempname,"(A3,I2.2,A8,I6.6,A4)") "Run", this%RunID,"_speck_t",tid,".out"
            fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
            call decomp_2d_write_plane(decompdir, tmpvar, dirid, jindx, fname, this%sp_gpC)

            if(this%isStratified) then
                jindx = jindx + 1 ! T
                write(tempname,"(A3,I2.2,A8,I6.6,A4)") "Run", this%RunID,"_specT_t",tid,".out"
                fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                call decomp_2d_write_plane(decompdir, tmpvar, dirid, jindx, fname, this%sp_gpC)
            endif

            if(this%fastCalcPressure .or. this%storePressure) then
                jindx = jindx + 1 ! p
                write(tempname,"(A3,I2.2,A8,I6.6,A4)") "Run", this%RunID,"_specp_t",tid,".out"
                fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                call decomp_2d_write_plane(decompdir, tmpvar, dirid, jindx, fname, this%sp_gpC)
            endif
        end if 

        nullify(S11_mean3D, S12_mean3D, S13_mean3D, S22_mean3D, S23_mean3D, S33_mean3D)

    end subroutine

    subroutine finalize_stats3D(this)
        class(igrid), intent(inout) :: this

        nullify(this%u_mean, this%v_mean, this%w_mean, this%uu_mean, this%uv_mean, this%uw_mean, this%vv_mean, this%vw_mean, this%ww_mean, this%disperuw_mean, this%dispervw_mean)
        nullify(this%mkeadv_mean, this%mkett_mean, this%mkedisp_mean, this%tkeadv_mean, this%tkett_mean, this%tkeprod_mean)
        nullify(this%u_mean3D, this%v_mean3D, this%w_mean3D, this%uu_mean3D, this%uv_mean3D, this%uw_mean3D, this%vv_mean3D, this%vw_mean3D, this%ww_mean3D, this%tketurbtranspx_mean3D, this%tketurbtranspy_mean3D, this%tketurbtranspz_mean3D)

        if(this%fastCalcPressure .or. this%storePressure) then
            nullify(this%p_mean,   this%mkept_mean,  this%tkept_mean)
            nullify(this%p_mean3D, this%pu_mean3D, this%pv_mean3D, this%pw_mean3D)
        endif

        if(.not. this%isInviscid) then
            nullify(this%mkevdif_mean, this%mkevdsp_mean, this%tkevdif_mean, this%tkevdsp_mean)
            nullify(this%viscdisp_mean3D, this%Siju1_mean3D, this%Siju2_mean3D, this%Siju3_mean3D)
        endif

        if(this%useSGS) then
           nullify(this%mkesgst_mean, this%mkesgsd_mean, this%tkesgst_mean, this%tkesgsd_mean)
           nullify(this%tau11_mean, this%tau12_mean, this%tau13_mean, this%tau22_mean, this%tau23_mean, this%tau33_mean)
           nullify(this%tau11_mean3D, this%tau12_mean3D, this%tau13_mean3D, this%tau22_mean3D, this%tau23_mean3D, this%tau33_mean3D)
           nullify(this%sgsdissp_mean3D, this%tauu1_mean3D, this%tauu2_mean3D, this%tauu3_mean3D)
        endif

        if(this%useWindTurbines) then
            nullify(this%turbfx_mean, this%turbfy_mean, this%turbfz_mean, this%tketurbf_mean, this%mketurbf_mean)
            nullify(this%turbfx_mean3D, this%turbfy_mean3D, this%turbfz_mean3D, this%uturbf_mean3D)
        endif

        if(this%isStratified) then
            nullify(this%T_mean, this%uT_mean, this%vT_mean, this%wT_mean, this%TT_mean)
            nullify(this%T_mean3D, this%uT_mean3D, this%vT_mean3D, this%wT_mean3D, this%TT_mean3D)
           if(this%useSGS) then
              nullify(this%q1_mean, this%q2_mean, this%q3_mean)
              nullify(this%q1_mean3D, this%q2_mean3D, this%q3_mean3D)
           endif
        endif

        deallocate(this%stats3D, this%horzavgstats, this%inst_horz_avg, this%runningSum_sc, this%xspectra_mean)
        if(this%useWindTurbines) deallocate(this%inst_horz_avg_turb, this%runningSum_sc_turb, this%runningSum_turb)
    end subroutine

    !--------------------------------Done 3D Statistics----------------------------------------------

    !!--------------------------------Beginning 1D Statistics----------------------------------------------
    !subroutine init_stats( this)
    !    class(igrid), intent(inout), target :: this
    !    type(decomp_info), pointer  :: gpC

    !    gpC => this%gpC
    !    this%tidSUM = 0

    !    if (this%isStratified) then
    !        allocate(this%zStats2dump(this%nz,33))
    !        allocate(this%runningSum(this%nz,33))
    !        allocate(this%TemporalMnNOW(this%nz,33))
    !        allocate(this%runningSum_sc(5))
    !        allocate(this%inst_horz_avg(5))
    !    else
    !        allocate(this%zStats2dump(this%nz,25))
    !        allocate(this%runningSum(this%nz,25))
    !        allocate(this%TemporalMnNOW(this%nz,25))
    !        allocate(this%runningSum_sc(3))
    !        allocate(this%inst_horz_avg(3))
    !    end if 

    !    if(this%useWindTurbines) then
    !        allocate(this%inst_horz_avg_turb(8*this%WindTurbineArr%nTurbines))
    !        allocate(this%runningSum_sc_turb(8*this%WindTurbineArr%nTurbines))
    !        allocate(this%runningSum_turb   (8*this%WindTurbineArr%nTurbines))
    !    endif

    !    ! mean velocities
    !    this%u_mean => this%zStats2dump(:,1);  this%v_mean  => this%zStats2dump(:,2);  this%w_mean => this%zStats2dump(:,3) 

    !    ! mean squared velocities
    !    this%uu_mean => this%zStats2dump(:,4); this%uv_mean => this%zStats2dump(:,5); this%uw_mean => this%zStats2dump(:,6)
    !                                           this%vv_mean => this%zStats2dump(:,7); this%vw_mean => this%zStats2dump(:,8) 
    !                                                                                  this%ww_mean => this%zStats2dump(:,9)

    !    ! SGS stresses
    !    this%tau11_mean => this%zStats2dump(:,10); this%tau12_mean => this%zStats2dump(:,11); this%tau13_mean => this%zStats2dump(:,12)
    !                                               this%tau22_mean => this%zStats2dump(:,13); this%tau23_mean => this%zStats2dump(:,14) 
    !                                                                                          this%tau33_mean => this%zStats2dump(:,15)

    !    ! SGS dissipation
    !    this%sgsdissp_mean => this%zStats2dump(:,16)

    !    ! velocity derivative products - for viscous dissipation
    !    this%viscdisp_mean => this%zStats2dump(:,17)

    !    ! means of velocity derivatives
    !    this%S11_mean => this%zStats2dump(:,18); this%S12_mean => this%zStats2dump(:,19); this%S13_mean => this%zStats2dump(:,20)
    !                                             this%S22_mean => this%zStats2dump(:,21); this%S23_mean => this%zStats2dump(:,22)
    !                                                                                      this%S33_mean => this%zStats2dump(:,23)

    !    ! SGS model coefficient
    !    this%sgscoeff_mean => this%zStats2dump(:,24)

    !    this%PhiM => this%zStats2dump(:,25)
    !    
    !    if (this%isStratified) then
    !        this%TT_mean => this%zStats2dump(:,30);  this%wT_mean => this%zStats2Dump(:,29);  this%vT_mean => this%zStats2Dump(:,28)
    !        this%uT_mean => this%zStats2dump(:,27);  this%T_mean => this%zStats2Dump(:,26); this%q1_mean => this%zStats2Dump(:,31)
    !        this%q2_mean => this%zStats2dump(:,32);  this%q3_mean => this%zStats2Dump(:,33)
    !    end if 
    !    this%runningSum_sc = zero
    !    this%runningSum = zero
    !    this%TemporalMnNOW = zero
    !    this%zStats2dump = zero
    !    this%inst_horz_avg = zero
    !    if(this%useWindTurbines) then
    !        this%inst_horz_avg_turb = zero
    !        this%runningSum_sc_turb = zero
    !        this%runningSum_turb    = zero
    !    endif
    !    nullify(gpC)
    !end subroutine

    !subroutine compute_stats(this)
    !    class(igrid), intent(inout), target :: this
    !    type(decomp_info), pointer :: gpC
    !    real(rkind), dimension(:,:,:), pointer :: rbuff1, rbuff2, rbuff3E, rbuff2E, rbuff3, rbuff4, rbuff5, rbuff5E, rbuff4E, rbuff6E, rbuff6!, rbuff7

    !    rbuff1  => this%rbuffxC(:,:,:,1); rbuff2  => this%rbuffyC(:,:,:,1);
    !    rbuff2E => this%rbuffyE(:,:,:,1); rbuff3E => this%rbuffzE(:,:,:,1);
    !    rbuff3 => this%rbuffzC(:,:,:,1); rbuff4E => this%rbuffzE(:,:,:,2);
    !    rbuff4 => this%rbuffzC(:,:,:,2); rbuff5E => this%rbuffzE(:,:,:,3)
    !    rbuff5 => this%rbuffzC(:,:,:,3); rbuff6E => this%rbuffzE(:,:,:,4)
    !    rbuff6 => this%rbuffzC(:,:,:,4); !rbuff7 => this%rbuffzC(:,:,:,5); 
    !    gpC => this%gpC

    !    this%tidSUM = this%tidSUM + 1


    !    ! Compute u - mean 
    !    call transpose_x_to_y(this%u,rbuff2,this%gpC)
    !    call transpose_y_to_z(rbuff2,rbuff3,this%gpC)
    !    call this%compute_z_mean(rbuff3, this%u_mean)
    !    !if (this%normByustar) this%u_mean = this%u_mean/this%sgsmodel%get_ustar()

    !    ! Compute v - mean 
    !    call transpose_x_to_y(this%v,rbuff2,this%gpC)
    !    call transpose_y_to_z(rbuff2,rbuff4,this%gpC)
    !    call this%compute_z_mean(rbuff4, this%v_mean)
    !    !if (this%normByustar)this%v_mean = this%v_mean/this%sgsmodel%get_ustar()

    !    ! Compute wC - mean 
    !    call transpose_x_to_y(this%wC,rbuff2,this%gpC)
    !    call transpose_y_to_z(rbuff2,rbuff5,this%gpC)
    !    call this%compute_z_mean(rbuff5, this%w_mean)
    !    !if (this%normByustar)this%w_mean = this%w_mean/this%sgsmodel%get_ustar()

    !    ! take w from x -> z decomp
    !    call transpose_x_to_y(this%w,rbuff2E,this%gpE)
    !    call transpose_y_to_z(rbuff2E,rbuff5E,this%gpE)

    !    ! take uE from x -> z decomp
    !    call transpose_x_to_y(this%uE,rbuff2E,this%gpE)
    !    call transpose_y_to_z(rbuff2E,rbuff3E,this%gpE)

    !    ! take vE from x -> z decomp
    !    call transpose_x_to_y(this%vE,rbuff2E,this%gpE)
    !    call transpose_y_to_z(rbuff2E,rbuff4E,this%gpE)

    !    ! uu mean
    !    rbuff6 = rbuff3*rbuff3
    !    call this%compute_z_mean(rbuff6, this%uu_mean)
    !    !if (this%normByustar)this%uu_mean = this%uu_mean/(this%sgsmodel%get_ustar()**2)

    !    ! uv mean
    !    rbuff6 = rbuff3*rbuff4
    !    call this%compute_z_mean(rbuff6, this%uv_mean)
    !    !if (this%normByustar)this%uv_mean = this%uv_mean/(this%sgsmodel%get_ustar()**2)

    !    ! uw mean
    !    rbuff6E = rbuff3E*rbuff5E
    !    call this%OpsPP%InterpZ_Edge2Cell(rbuff6E,rbuff6)
    !    call this%compute_z_mean(rbuff6, this%uw_mean)
    !    !if (this%normByustar)this%uw_mean = this%uw_mean/(this%sgsmodel%get_ustar()**2)

    !    ! vv mean 
    !    rbuff6 = rbuff4*rbuff4
    !    call this%compute_z_mean(rbuff6, this%vv_mean)
    !    !if (this%normByustar)this%vv_mean = this%vv_mean/(this%sgsmodel%get_ustar()**2)

    !    ! vw mean 
    !    rbuff6E = rbuff4E*rbuff5E
    !    call this%OpsPP%InterpZ_Edge2Cell(rbuff6E,rbuff6)
    !    call this%compute_z_mean(rbuff6, this%vw_mean)
    !    !if (this%normByustar)this%vw_mean = this%vw_mean/(this%sgsmodel%get_ustar()**2)

    !    ! ww mean 
    !    rbuff6 = rbuff5*rbuff5
    !    call this%compute_z_mean(rbuff6, this%ww_mean)
    !    !if (this%normByustar)this%ww_mean = this%ww_mean/(this%sgsmodel%get_ustar()**2)

    !    ! Statified Stuff
    !    if (this%isStratified) then
    !        ! T mean
    !        call transpose_x_to_y(this%T,rbuff2,this%gpC)
    !        call transpose_y_to_z(rbuff2,rbuff6,this%gpC)
    !        call this%compute_z_mean(rbuff6, this%T_mean)

    !        ! uT mean
    !        rbuff3 = rbuff3*rbuff6
    !        call this%compute_z_mean(rbuff3, this%uT_mean)
    !        
    !        ! vT mean
    !        rbuff4 = rbuff4*rbuff6
    !        call this%compute_z_mean(rbuff4, this%vT_mean)
    !        
    !        ! wT mean
    !        rbuff5 = rbuff5*rbuff6
    !        call this%compute_z_mean(rbuff5, this%wT_mean)
    !        
    !        ! TT mean
    !        rbuff6 = rbuff6*rbuff6
    !        call this%compute_z_mean(rbuff6, this%TT_mean)
    !    end if 


    !    if (this%useSGS) then
    !        ! tau_11
    !        call transpose_x_to_y(this%tauSGS_ij(:,:,:,1),rbuff2,this%gpC)
    !        call transpose_y_to_z(rbuff2,rbuff3,this%gpC)
    !        call this%compute_z_mean(rbuff3, this%tau11_mean)
    !        if (this%normByustar)this%tau11_mean = this%tau11_mean/(this%sgsmodel%get_ustar()**2)

    !        ! tau_12
    !        call transpose_x_to_y(this%tauSGS_ij(:,:,:,2),rbuff2,this%gpC)
    !        call transpose_y_to_z(rbuff2,rbuff3,this%gpC)
    !        call this%compute_z_mean(rbuff3, this%tau12_mean)
    !        if (this%normByustar)this%tau12_mean = this%tau12_mean/(this%sgsmodel%get_ustar()**2)

    !        ! tau_13
    !        call transpose_x_to_y(this%tau13,rbuff2E,this%gpE)
    !        call transpose_y_to_z(rbuff2E,rbuff3E,this%gpE)
    !        rbuff3E(:,:,1) = -(this%sgsmodel%get_ustar()**2)
    !        call this%OpsPP%InterpZ_Edge2Cell(rbuff3E,rbuff3)
    !        call transpose_z_to_y(rbuff3,rbuff2,this%gpC)
    !        call transpose_y_to_x(rbuff2,this%tauSGS_ij(:,:,:,3),this%gpC)
    !        call this%compute_z_mean(rbuff3, this%tau13_mean)
    !        if (this%normByustar)this%tau13_mean = this%tau13_mean/(this%sgsmodel%get_ustar()**2)

    !        ! tau_22
    !        call transpose_x_to_y(this%tauSGS_ij(:,:,:,4),rbuff2,this%gpC)
    !        call transpose_y_to_z(rbuff2,rbuff3,this%gpC)
    !        call this%compute_z_mean(rbuff3, this%tau22_mean)
    !        if (this%normByustar)this%tau22_mean = this%tau22_mean/(this%sgsmodel%get_ustar()**2)

    !        ! tau_23
    !        call transpose_x_to_y(this%tau23,rbuff2E,this%gpE)
    !        call transpose_y_to_z(rbuff2E,rbuff3E,this%gpE)
    !        rbuff3E(:,:,1) = -(this%sgsmodel%get_ustar()**2)*this%sgsmodel%get_vmean()/this%sgsmodel%get_umean()
    !        call this%OpsPP%InterpZ_Edge2Cell(rbuff3E,rbuff3)
    !        call transpose_z_to_y(rbuff3,rbuff2,this%gpC)
    !        call transpose_y_to_x(rbuff2,this%tauSGS_ij(:,:,:,5),this%gpC)
    !        call this%compute_z_mean(rbuff3, this%tau23_mean)
    !        if (this%normByustar)this%tau23_mean = this%tau23_mean/(this%sgsmodel%get_ustar()**2)

    !        ! tau_33
    !        call transpose_x_to_y(this%tauSGS_ij(:,:,:,6),rbuff2,this%gpC)
    !        call transpose_y_to_z(rbuff2,rbuff3,this%gpC)
    !        call this%compute_z_mean(rbuff3, this%tau33_mean)
    !        if (this%normByustar)this%tau33_mean = this%tau33_mean/(this%sgsmodel%get_ustar()**2)


    !        ! sgs dissipation
    !        rbuff1 = this%tauSGS_ij(:,:,:,1)*this%tauSGS_ij(:,:,:,1) + &
    !                 this%tauSGS_ij(:,:,:,2)*this%tauSGS_ij(:,:,:,2) + &
    !                 this%tauSGS_ij(:,:,:,3)*this%tauSGS_ij(:,:,:,3)
    !        rbuff1 = rbuff1 + two*(this%tauSGS_ij(:,:,:,4)*this%tauSGS_ij(:,:,:,4) + &
    !                               this%tauSGS_ij(:,:,:,5)*this%tauSGS_ij(:,:,:,5) + &
    !                               this%tauSGS_ij(:,:,:,6)*this%tauSGS_ij(:,:,:,6) )
    !        rbuff1 = rbuff1/(this%nu_SGS + 1.0d-14)         ! note: factor of half is in dump_stats

    !        call transpose_x_to_y(rbuff1,rbuff2,this%gpC)
    !        call transpose_y_to_z(rbuff2,rbuff3,this%gpC)
    !        call this%compute_z_mean(rbuff3, this%sgsdissp_mean)

    !        ! viscous dissipation- *****????? Is rbuff1 contaminated after transpose_x_to_y? *****?????
    !        rbuff1 = rbuff1/(this%nu_SGS + 1.0d-14)        ! note: factor of fourth is in dump_stats
    !        call transpose_x_to_y(rbuff1,rbuff2,this%gpC)
    !        call transpose_y_to_z(rbuff2,rbuff3,this%gpC)
    !        call this%compute_z_mean(rbuff3, this%viscdisp_mean)

    !        ! note: factor of half in all S_** is in dump_stats
    !        ! S_11
    !        rbuff1 = this%tauSGS_ij(:,:,:,1)/(this%nu_SGS + 1.0d-14)
    !        call transpose_x_to_y(rbuff1,rbuff2,this%gpC)
    !        call transpose_y_to_z(rbuff2,rbuff3,this%gpC)
    !        call this%compute_z_mean(rbuff3, this%S11_mean)

    !        ! S_12
    !        rbuff1 = this%tauSGS_ij(:,:,:,2)/(this%nu_SGS + 1.0d-14)
    !        call transpose_x_to_y(rbuff1,rbuff2,this%gpC)
    !        call transpose_y_to_z(rbuff2,rbuff3,this%gpC)
    !        call this%compute_z_mean(rbuff3, this%S12_mean)

    !        ! S_13
    !        rbuff1 = this%tauSGS_ij(:,:,:,3)/(this%nu_SGS + 1.0d-14)
    !        call transpose_x_to_y(rbuff1,rbuff2,this%gpC)
    !        call transpose_y_to_z(rbuff2,rbuff3,this%gpC)
    !        call this%compute_z_mean(rbuff3, this%S13_mean)

    !        ! S_22
    !        rbuff1 = this%tauSGS_ij(:,:,:,4)/(this%nu_SGS + 1.0d-14)
    !        call transpose_x_to_y(rbuff1,rbuff2,this%gpC)
    !        call transpose_y_to_z(rbuff2,rbuff3,this%gpC)
    !        call this%compute_z_mean(rbuff3, this%S22_mean)

    !        ! S_23
    !        rbuff1 = this%tauSGS_ij(:,:,:,5)/(this%nu_SGS + 1.0d-14)
    !        call transpose_x_to_y(rbuff1,rbuff2,this%gpC)
    !        call transpose_y_to_z(rbuff2,rbuff3,this%gpC)
    !        call this%compute_z_mean(rbuff3, this%S23_mean)

    !        ! S_33
    !        rbuff1 = this%tauSGS_ij(:,:,:,6)/(this%nu_SGS + 1.0d-14)
    !        call transpose_x_to_y(rbuff1,rbuff2,this%gpC)
    !        call transpose_y_to_z(rbuff2,rbuff3,this%gpC)
    !        call this%compute_z_mean(rbuff3, this%S33_mean)

    !        ! sgs coefficient
    !        call transpose_x_to_y(this%c_SGS,rbuff2,this%gpC)
    !        call transpose_y_to_z(rbuff2,rbuff3,this%gpC)
    !        !call this%compute_z_mean(rbuff3, this%sgscoeff_mean)    ! -- averaging not needed
    !        this%sgscoeff_mean(:) = rbuff3(1,1,:)
    !   
    !        
    !        if (this%isStratified) then
    !            ! q1
    !            call transpose_x_to_y(this%q1,rbuff2,this%gpC)
    !            call transpose_y_to_z(rbuff2,rbuff3,this%gpC)
    !            call this%compute_z_mean(rbuff3, this%q1_mean)    

    !            ! q2
    !            call transpose_x_to_y(this%q2,rbuff2,this%gpC)
    !            call transpose_y_to_z(rbuff2,rbuff3,this%gpC)
    !            call this%compute_z_mean(rbuff3, this%q2_mean)    

    !            ! q3
    !            call transpose_x_to_y(this%q3,rbuff2E,this%gpE)
    !            call transpose_y_to_z(rbuff2E,rbuff3E,this%gpE)
    !            rbuff3E(:,:,1) = this%wTh_surf
    !            call this%OpsPP%InterpZ_Edge2Cell(rbuff3E,rbuff3)
    !            call this%compute_z_mean(rbuff3, this%q3_mean)    
    !        end if 
    !    end if

    !    rbuff1 = this%duidxjC(:,:,:,3)*this%mesh(:,:,:,3)
    !    call transpose_x_to_y(rbuff1,rbuff2,this%gpC)
    !    call transpose_y_to_z(rbuff2,rbuff3,this%gpC)
    !    call this%compute_z_mean(rbuff3, this%PhiM)
    !    if (this%useSGS) then
    !     this%PhiM = this%PhiM*kappa/this%sgsmodel%get_ustar()
    !    end if

    !    this%runningSum = this%runningSum + this%zStats2dump

    !    !write(*,*) 'In stats'
    !    !write(*,*) 'umean', maxval(this%u_mean), minval(this%u_mean)
    !    !write(*,*) 'vmean', maxval(this%v_mean), minval(this%v_mean)
    !    !write(*,*) 'wmean', maxval(this%w_mean), minval(this%w_mean)

    !    ! instantaneous horizontal averages of some quantities
    !    if (this%useSGS) this%inst_horz_avg(1) = this%sgsmodel%get_ustar()
    !    ! this%inst_horz(2) and (3) are computed in getRHS_SGS_WallM
    !    if(this%isStratified) then
    !        this%inst_horz_avg(4) = this%invObLength
    !        this%inst_horz_avg(5) = this%wTh_surf
    !    endif
    !    ! this%inst_horz_avg_turb(1:5*this%WindTurbineArr%nTurbines) is computed in this%WindTurbineArr%getForceRHS
    !    this%runningSum_sc = this%runningSum_sc + this%inst_horz_avg
    !    if(this%useWindTurbines) this%runningSum_sc_turb = this%runningSum_sc_turb + this%inst_horz_avg_turb

    !end subroutine 

    !subroutine dump_stats(this)
    !    use basic_io, only: write_2d_ascii, write_2D_binary
    !    use exits, only: message
    !    use kind_parameters, only: clen, mpirkind
    !    use mpi
    !    class(igrid), intent(inout), target :: this
    !    character(len=clen) :: fname
    !    character(len=clen) :: tempname
    !    integer :: tid, ierr

    !    this%TemporalMnNOW = this%runningSum/real(this%tidSUM,rkind)
    !    tid = this%step

    !    ! compute (u_i'u_j')
    !    this%TemporalMnNOW(:,4) = this%TemporalMnNOW(:,4) - this%TemporalMnNOW(:,1)*this%TemporalMnNOW(:,1)
    !    this%TemporalMnNOW(:,5) = this%TemporalMnNOW(:,5) - this%TemporalMnNOW(:,1)*this%TemporalMnNOW(:,2)
    !    this%TemporalMnNOW(:,6) = this%TemporalMnNOW(:,6) - this%TemporalMnNOW(:,1)*this%TemporalMnNOW(:,3)
    !    this%TemporalMnNOW(:,7) = this%TemporalMnNOW(:,7) - this%TemporalMnNOW(:,2)*this%TemporalMnNOW(:,2)
    !    this%TemporalMnNOW(:,8) = this%TemporalMnNOW(:,8) - this%TemporalMnNOW(:,2)*this%TemporalMnNOW(:,3)
    !    this%TemporalMnNOW(:,9) = this%TemporalMnNOW(:,9) - this%TemporalMnNOW(:,3)*this%TemporalMnNOW(:,3)

    !    if (this%isStratified) then
    !        this%TemporalMnNOW(:,27) = this%TemporalMnNOW(:,27) - this%TemporalMnNOW(:,1)*this%TemporalMnNOW(:,26)
    !        this%TemporalMnNOW(:,28) = this%TemporalMnNOW(:,28) - this%TemporalMnNOW(:,2)*this%TemporalMnNOW(:,26)
    !        this%TemporalMnNOW(:,29) = this%TemporalMnNOW(:,29) - this%TemporalMnNOW(:,3)*this%TemporalMnNOW(:,26)
    !        this%TemporalMnNOW(:,30) = this%TemporalMnNOW(:,30) - this%TemporalMnNOW(:,26)*this%TemporalMnNOW(:,26)
    !    end if 

    !    ! compute sgs dissipation
    !    this%TemporalMnNOW(:,16) = half*this%TemporalMnNOW(:,16)

    !    ! compute viscous dissipation
    !    this%TemporalMnNOW(:,17) = this%TemporalMnNOW(:,17) - (                        &
    !                               this%TemporalMnNOW(:,18)*this%TemporalMnNOW(:,18) + &
    !                               this%TemporalMnNOW(:,21)*this%TemporalMnNOW(:,21) + &
    !                               this%TemporalMnNOW(:,23)*this%TemporalMnNOW(:,23) + &
    !                          two*(this%TemporalMnNOW(:,19)*this%TemporalMnNOW(:,19) + &
    !                               this%TemporalMnNOW(:,20)*this%TemporalMnNOW(:,20) + & 
    !                               this%TemporalMnNOW(:,22)*this%TemporalMnNOW(:,22)))
    !    this%TemporalMnNOW(:,17) = half*this%TemporalMnNOW(:,17)/this%Re     ! note: this is actually 2/Re*(..)/4

    !    if(this%useWindTurbines) then
    !        this%runningSum_turb = zero
    !        call MPI_reduce(this%runningSum_sc_turb, this%runningSum_turb, 8*this%WindTurbineArr%nTurbines, mpirkind, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
    !    endif
    !    if (nrank == 0) then
    !        write(tempname,"(A3,I2.2,A2,I6.6,A4)") "Run", this%RunID,"_t",tid,".stt"
    !        fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
    !        call write_2d_ascii(this%TemporalMnNOW,fname)
    !        !call write_2D_binary(TemporalMnNOW,fname)

    !        write(tempname,"(A3,I2.2,A2,I6.6,A4)") "Run", this%RunID,"_t",tid,".sth"   ! time and horz averages of scalars
    !        fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
    !        open(unit=771,file=fname,status='unknown')
    !        if(this%useWindTurbines) then
    !            write(771,'(e19.12,1x,i7,1x,8008(e19.12,1x))') this%tsim, this%tidSUM, this%runningSum_sc/real(this%tidSUM,rkind), this%runningSum_turb/real(this%tidSUM,rkind) ! change if using more than 1000 turbines
    !        else
    !            write(771,'(e19.12,1x,i7,1x,5(e19.12,1x))') this%tsim, this%tidSUM, this%runningSum_sc/real(this%tidSUM,rkind)
    !        endif
    !        close(771)
    !    end if
    !    call message(1, "Just dumped a .stt file")
    !    call message(2, "Number ot tsteps averaged:",this%tidSUM)

    !end subroutine

    !subroutine compute_z_fluct(this,fin)
    !    use reductions, only: P_SUM
    !    class(igrid), intent(in), target :: this
    !    real(rkind), dimension(:,:,:), intent(inout) :: fin
    !    integer :: k
    !    real(rkind) :: fmean

    !    do k = 1,size(fin,3)
    !        fmean = P_SUM(sum(fin(:,:,k)))/(real(this%nx,rkind)*real(this%ny,rkind))
    !        fin(:,:,k) = fin(:,:,k) - fmean
    !    end do 

    !end subroutine

    !subroutine compute_y_mean(this, arr_in, arr_out)
    !    use reductions, only: P_SUM
    !    class(igrid), intent(in), target :: this
    !    real(rkind), dimension(:,:,:), intent(in) :: arr_in
    !    real(rkind), dimension(:,:), intent(out) :: arr_out
    !    integer :: k, i

    !    do k = 1,size(arr_in,3)
    !      do i = 1,size(arr_in,1)
    !        arr_out(i,k) = P_SUM(sum(arr_in(:,:,k)))/(real(this%nx,rkind)*real(this%ny,rkind))
    !    end do 

    !end subroutine

    subroutine compute_z_mean(this, arr_in, vec_out)
        use reductions, only: P_SUM
        class(igrid), intent(in), target :: this
        real(rkind), dimension(:,:,:), intent(in) :: arr_in
        real(rkind), dimension(:), intent(out) :: vec_out
        integer :: k

        do k = 1,size(arr_in,3)
            vec_out(k) = P_SUM(sum(arr_in(:,:,k)))/(real(this%nx,rkind)*real(this%ny,rkind))
        end do 

    end subroutine

    !subroutine finalize_stats(this)
    !    class(igrid), intent(inout) :: this

    !    nullify(this%u_mean, this%v_mean, this%w_mean, this%uu_mean, this%uv_mean, this%uw_mean, this%vv_mean, this%vw_mean, this%ww_mean)
    !    nullify(this%tau11_mean, this%tau12_mean, this%tau13_mean, this%tau22_mean, this%tau23_mean, this%tau33_mean)
    !    nullify(this%S11_mean, this%S12_mean, this%S13_mean, this%S22_mean, this%S23_mean, this%S33_mean)
    !    nullify(this%sgsdissp_mean, this%viscdisp_mean, this%sgscoeff_mean)
    !    if (allocated(this%zStats2dump)) deallocate(this%zStats2dump)
    !    if (allocated(this%runningSum)) deallocate(this%runningSum)
    !    if (allocated(this%TemporalMnNow)) deallocate(this%TemporalMnNOW)
    !    if (allocated(this%runningSum_sc)) deallocate(this%runningSum_sc)
    !    if(this%useWindTurbines) deallocate(this%inst_horz_avg_turb, this%runningSum_sc_turb, this%runningSum_turb)
    !end subroutine 

    !!--------------------------------Done 1D Statistics----------------------------------------------

    subroutine dump_pointProbes(this)
        use kind_parameters, only: mpirkind
        !class(igrid), intent(inout) :: this
        class(igrid), intent(in) :: this
        !character(len=clen) :: fname
        !character(len=clen) :: tempname
        !integer :: ierr

        ! ADITYA -> NIRANJAN: Why isn't this quantity inside turbarray? It
        ! segfaults for certain specific casses. 
        if(this%useWindTurbines) then
         !   this%runningSum_turb = zero
            !call MPI_reduce(this%inst_horz_avg_turb, this%runningSum_turb, 8*this%WindTurbineArr%nTurbines, mpirkind, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
            !if(nrank == 0) then
            !    write(tempname,"(A3,I2.2,A15)") "Run", this%RunID,"_timeseries.prb"
            !    fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
            !    open(unit=10,file=fname,status='old',action='write',position='append',iostat=ierr)
            !    if(ierr .ne. 0) open(unit=10,file=fname,status='replace')
            !    write(10,'(1000(e19.12,1x))') this%tsim, this%inst_horz_avg, this%runningSum_turb
            !    close(10)
            !end if
        endif

    end subroutine 

    subroutine dump_planes(this)
        use decomp_2d_io
        class(igrid), intent(in) :: this
        integer :: nxplanes, nyplanes, nzplanes
        integer :: idx, pid, dirid, tid
        character(len=clen) :: fname
        character(len=clen) :: tempname



        tid = this%step 
        if (allocated(this%xplanes)) then
            nxplanes = size(this%xplanes)
            dirid = 1
            do idx = 1,nxplanes
                pid = this%xplanes(idx)
            
                write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_x",pid,".plu"
                fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                call decomp_2d_write_plane(1,this%u,dirid, pid, fname, this%gpC)

                write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_x",pid,".plv"
                fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                call decomp_2d_write_plane(1,this%v,dirid, pid, fname, this%gpC)

                write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_x",pid,".plw"
                fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                call decomp_2d_write_plane(1,this%wC,dirid, pid, fname, this%gpC)
                
                if (this%isStratified) then
                    write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_x",pid,".plT"
                    fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                    call decomp_2d_write_plane(1,this%T,dirid, pid, fname, this%gpC)
                end if

                if (this%fastCalcPressure) then
                    write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_x",pid,".plP"
                    fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                    call decomp_2d_write_plane(1,this%Pressure,dirid, pid, fname, this%gpC)
                end if 


                if (this%computevorticity) then
                    write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_x",pid,".pox"
                    fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                    call decomp_2d_write_plane(1,this%ox,dirid, pid, fname, this%gpC)
                    
                    write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_x",pid,".poy"
                    fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                    call decomp_2d_write_plane(1,this%oy,dirid, pid, fname, this%gpC)
                    
                    write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_x",pid,".poz"
                    fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                    call decomp_2d_write_plane(1,this%oz,dirid, pid, fname, this%gpC)
                end if 

                ! planes for KS preprocess
                if (this%PreProcessForKS) then
                    write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_x",pid,".ksu"
                    fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                    call decomp_2d_write_plane(1,this%uFil4KS,dirid, pid, fname, this%gpC)

                    write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_x",pid,".ksv"
                    fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                    call decomp_2d_write_plane(1,this%vFil4KS,dirid, pid, fname, this%gpC)

                    write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_x",pid,".ksw"
                    fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                    call decomp_2d_write_plane(1,this%wFil4KS,dirid, pid, fname, this%gpC)
                end if 
            end do 
        end if 
            
            
        if (allocated(this%yplanes)) then
            nyplanes = size(this%yplanes)
            dirid = 2
            do idx = 1,nyplanes
                pid = this%yplanes(idx)

                write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_y",pid,".plu"
                fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                call decomp_2d_write_plane(1,this%u,dirid, pid, fname, this%gpC)

                write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_y",pid,".plv"
                fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                call decomp_2d_write_plane(1,this%v,dirid, pid, fname, this%gpC)

                write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_y",pid,".plw"
                fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                call decomp_2d_write_plane(1,this%wC,dirid, pid, fname, this%gpC)
                
                if (this%isStratified) then
                    write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_y",pid,".plT"
                    fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                    call decomp_2d_write_plane(1,this%T,dirid, pid, fname, this%gpC)
                end if
                
                if (this%fastCalcPressure) then
                    write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_y",pid,".plP"
                    fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                    call decomp_2d_write_plane(1,this%Pressure,dirid, pid, fname, this%gpC)
                end if 

                if (this%computevorticity) then
                    write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_y",pid,".pox"
                    fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                    call decomp_2d_write_plane(1,this%ox,dirid, pid, fname, this%gpC)
                    
                    write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_y",pid,".poy"
                    fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                    call decomp_2d_write_plane(1,this%oy,dirid, pid, fname, this%gpC)
                    
                    write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_y",pid,".poz"
                    fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                    call decomp_2d_write_plane(1,this%oz,dirid, pid, fname, this%gpC)
                end if 
                
                ! planes for KS preprocess
                if (this%PreProcessForKS) then
                    write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_y",pid,".ksu"
                    fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                    call decomp_2d_write_plane(1,this%uFil4KS,dirid, pid, fname, this%gpC)

                    write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_y",pid,".ksv"
                    fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                    call decomp_2d_write_plane(1,this%vFil4KS,dirid, pid, fname, this%gpC)

                    write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_y",pid,".ksw"
                    fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                    call decomp_2d_write_plane(1,this%wFil4KS,dirid, pid, fname, this%gpC)

                end if

            end do 
        end if 
        
        
        if (allocated(this%zplanes)) then
            nzplanes = size(this%zplanes)
            dirid = 3
            do idx = 1,nzplanes
                pid = this%zplanes(idx)

                write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_z",pid,".plu"
                fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                call decomp_2d_write_plane(1,this%u,dirid, pid, fname, this%gpC)

                write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_z",pid,".plv"
                fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                call decomp_2d_write_plane(1,this%v,dirid, pid, fname, this%gpC)

                write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_z",pid,".plw"
                fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                call decomp_2d_write_plane(1,this%wC,dirid, pid, fname, this%gpC)
                
                if (this%isStratified) then
                    write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_z",pid,".plT"
                    fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                    call decomp_2d_write_plane(1,this%T,dirid, pid, fname, this%gpC)
                end if
                
                if (this%fastCalcPressure) then
                    write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_z",pid,".plP"
                    fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                    call decomp_2d_write_plane(1,this%Pressure,dirid, pid, fname, this%gpC)
                end if 

                if (this%computevorticity) then
                    write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_x",pid,".pox"
                    fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                    call decomp_2d_write_plane(1,this%ox,dirid, pid, fname, this%gpC)
                    
                    write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_x",pid,".poy"
                    fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                    call decomp_2d_write_plane(1,this%oy,dirid, pid, fname, this%gpC)
                    
                    write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_x",pid,".poz"
                    fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                    call decomp_2d_write_plane(1,this%oz,dirid, pid, fname, this%gpC)
                end if 

                ! planes for KS preprocess
                if (this%PreProcessForKS) then
                    write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_z",pid,".ksu"
                    fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                    call decomp_2d_write_plane(1,this%uFil4KS,dirid, pid, fname, this%gpC)

                    write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_z",pid,".ksv"
                    fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                    call decomp_2d_write_plane(1,this%vFil4KS,dirid, pid, fname, this%gpC)

                    write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_z",pid,".ksw"
                    fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                    call decomp_2d_write_plane(1,this%wFil4KS,dirid, pid, fname, this%gpC)
                end if
            end do 
        end if 
        call message(1, "Dumped Planes.")        
    end subroutine 



    !subroutine getfilteredSpeedSqAtWall(this)
    !    class(igrid), intent(inout), target :: this

    !    real(rkind), dimension(:,:,:), pointer :: rbuffx1, rbuffx2
    !    complex(rkind), dimension(:,:,:), pointer :: cbuffy, tauWallH

    !    if (this%useSGS) then
    !        cbuffy => this%cbuffyC(:,:,:,1); tauWallH => this%cbuffzC(:,:,:,1)     
    !        rbuffx1 => this%filteredSpeedSq; rbuffx2 => this%rbuffxC(:,:,:,1)

    !        call transpose_y_to_z(this%uhat,tauWallH,this%sp_gpC)
    !        call this%spectC%SurfaceFilter_ip(tauWallH(:,:,1))
    !        call transpose_z_to_y(tauWallH,cbuffy, this%sp_gpC)
    !        call this%spectC%ifft(cbuffy,rbuffx1)

    !        call transpose_y_to_z(this%vhat,tauWallH,this%sp_gpC)
    !        call this%spectC%SurfaceFilter_ip(tauWallH(:,:,1))
    !        call transpose_z_to_y(tauWallH,cbuffy, this%sp_gpC)
    !        call this%spectC%ifft(cbuffy,rbuffx2)

    !        rbuffx1 = rbuffx1*rbuffx1
    !        rbuffx2 = rbuffx2*rbuffx2
    !        rbuffx1 = rbuffx1 + rbuffx2
    !    end if 

    !end subroutine  


    subroutine start_io(this, dumpInitField)
        class(igrid), target, intent(inout) :: this 
        character(len=clen) :: fname
        character(len=clen) :: tempname
        !character(len=clen) :: command
        character(len=clen) :: OutputDir
        !integer :: system 
        integer :: runIDX 
        logical :: isThere
        integer :: tag, idx, status(MPI_STATUS_SIZE), ierr
        integer, dimension(:,:), allocatable        :: xst,xen,xsz
        logical, optional, intent(in) :: dumpInitField

        ! Create data sharing info
        !if (nrank == 0) then
            allocate(xst(0:nproc-1,3),xen(0:nproc-1,3),xsz(0:nproc-1,3))
            xst = 0; xen = 0; xsz = 0;
        !end if


        ! communicate local processor grid info (Assume x-decomposition)
        if (nrank == 0) then
            xst(0,:) = this%gpC%xst
            xen(0,:) = this%gpC%xen
            
            tag = 0
            do idx = 1,nproc-1
                call MPI_RECV(xst(idx,:), 3, MPI_INTEGER, idx, tag,&
                    MPI_COMM_WORLD, status, ierr)
            end do 
            tag = 1
            do idx = 1,nproc-1
                call MPI_RECV(xen(idx,:), 3, MPI_INTEGER, idx, tag,&
                    MPI_COMM_WORLD, status, ierr)
            end do
           tag = 2
            do idx = 1,nproc-1
                call MPI_RECV(xsz(idx,:), 3, MPI_INTEGER, idx, tag,&
                    MPI_COMM_WORLD, status, ierr)
            end do

        else
            tag = 0
            call MPI_SEND(this%gpC%xst, 3, MPI_INTEGER, 0, tag, &
                 &      MPI_COMM_WORLD, ierr)
            tag = 1
            call MPI_SEND(this%gpC%xen, 3, MPI_INTEGER, 0, tag, &
                 &      MPI_COMM_WORLD, ierr)
            tag = 2
            call MPI_SEND(this%gpC%xsz, 3, MPI_INTEGER, 0, tag, &
                 &      MPI_COMM_WORLD, ierr)

        end if 

        OutputDir = this%outputdir
        runIDX = this%runID
        
        inquire(FILE=trim(OutputDir), exist=isThere)
        if (nrank == 0) then
            write(tempname,"(A3,I2.2,A6,A4)") "Run", runIDX, "HEADER",".txt"
            fname = OutputDir(:len_trim(OutputDir))//"/"//trim(tempname)

            open (this%headerfid, file=trim(fname), FORM='formatted', STATUS='replace',ACTION='write')
            write(this%headerfid,*)"========================================================================="
            write(this%headerfid,*)"---------------------  Header file for MATLAB ---------------------------"
            write(this%headerfid,"(A9,A10,A10,A10,A10,A10,A10)") "PROC", "xst", "xen", "yst", "yen","zst","zen"
            write(this%headerfid,*)"-------------------------------------------------------------------------"
            do idx = 0,nproc-1
                write(this%headerfid,"(I8,6I10)") idx, xst(idx,1), xen(idx,1), xst(idx,2), xen(idx,2), xst(idx,3), xen(idx,3)
            end do 
            write(this%headerfid,*)"-------------------------------------------------------------------------"
            write(this%headerfid,*)"Dumps made at:"
        end if
        call mpi_barrier(mpi_comm_world,ierr)
        
        !if (nrank == 0) then
            deallocate(xst, xen, xsz)
        !end if 

        if (present(dumpInitField)) then
            if (dumpInitField) then
                call message(0,"Performing initialization data dump.")
                call this%dumpFullField(this%u,'uVel')
                call this%dumpFullField(this%v,'vVel')
                call this%dumpFullField(this%wC,'wVel')
                call this%dumpVisualizationInfo()
                if (this%isStratified .or. this%initspinup) call this%dumpFullField(this%T,'potT')
                if (this%fastCalcPressure) call this%dumpFullField(this%pressure,'prss')
                if (this%useWindTurbines) then
                    this%WindTurbineArr%dumpTurbField = .true.
                    this%WindTurbineArr%step = this%step-1
                endif
                call message(0,"Done with the initialization data dump.")
            end if
        end if
    end subroutine
    
    subroutine finalize_io(this)
        class(igrid), intent(in) :: this

        if (nrank == 0) then
            write(this%headerfid,*) "--------------------------------------------------------------"
            write(this%headerfid,*) "------------------ END OF HEADER FILE ------------------------"
            close(this%headerfid)
        end if 
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
            uBC_bottom      = -1; vBC_bottom      = -1;
            dUdzBC_bottom   =  0; dVdzBC_bottom   =  0;
            WdUdzBC_bottom  =  0; WdVdzBC_bottom  =  0;
            UWBC_bottom     = +1; VWBC_bottom     = +1;
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
            uBC_top      = -1; vBC_top      = -1;
            dUdzBC_top   =  0; dVdzBC_top   =  0;
            WdUdzBC_top  =  0; WdVdzBC_top  =  0;
            UWBC_top     = +1; VWBC_top     = +1;   
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
         end select 
         select case (botBC_Temp)
         case (0) ! Dirichlet BC for temperature at the bottom
            TBC_bottom = 0; dTdzBC_bottom = 0; WTBC_bottom = -1; 
            WdTdzBC_bottom = 0;      
         case(1)  ! Homogenenous Neumann BC at the bottom
            TBC_bottom = 1; dTdzBC_bottom = -1; WTBC_bottom = -1;
            WdTdzBC_bottom = 1
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

    subroutine correctPressureRotationalForm(this)
        class(igrid), intent(inout) :: this 
        real(rkind) :: mfact, meanK

        this%rbuffxC(:,:,:,1) = 0.5d0*(this%u*this%u + this%v*this%v + this%wC*this%wC)
        mfact = one/(real(this%nx,rkind)*real(this%ny,rkind)*real(this%nz,rkind))
        meanK = p_sum(sum(this%rbuffxC(:,:,:,1)))*mfact
        this%pressure = this%pressure + meanK
        this%pressure = this%pressure - this%rbuffxC(:,:,:,1)
    end subroutine

    function get_dt(this, recompute) result(val)
        class(igrid), intent(inout) :: this 
        logical, intent(in), optional :: recompute
        real(rkind) :: val

        if (present(recompute)) then
           if (recompute) call this%compute_deltaT
        end if
        val = this%dt
    end function

    subroutine interpolate_cellField_to_edgeField(this, rxC, rxE, bc1, bc2)
        class(igrid), intent(inout) :: this 
        real(rkind), intent(in),  dimension(this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3)) :: rxC 
        real(rkind), intent(out), dimension(this%gpE%xsz(1), this%gpE%xsz(2), this%gpE%xsz(3)) :: rxE
        integer, intent(in) :: bc1, bc2

        call transpose_x_to_y(rxC, this%rbuffyC(:,:,:,1), this%gpC)
        call transpose_y_to_z(this%rbuffyC(:,:,:,1), this%rbuffzC(:,:,:,1), this%gpC)
        call this%Pade6opZ%interpz_E2C(this%rbuffzC(:,:,:,1), this%rbuffzE(:,:,:,1), bc1,bc2) 
        call transpose_z_to_y(this%rbuffzE(:,:,:,1), this%rbuffyE(:,:,:,1), this%gpE)
        call transpose_y_to_x(this%rbuffyE(:,:,:,1), rxE, this%gpE)

    end subroutine 
end module 
