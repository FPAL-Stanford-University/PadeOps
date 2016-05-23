module IncompressibleGridWallM
    use kind_parameters, only: rkind, clen
    use constants, only: imi, zero,one,two,three,half,fourth, pi, kappa 
    use GridMod, only: grid
    use gridtools, only: alloc_buffs, destroy_buffs
    use igrid_hooks, only: meshgen_WallM, initfields_wallM, set_planes_io
    use decomp_2d
    use StaggOpsMod, only: staggOps  
    use exits, only: GracefulExit, message
    use spectralMod, only: spectral  
    use PoissonMod, only: poisson
    use mpi 
    use reductions, only: p_maxval, p_sum
    use timer, only: tic, toc
    use PadePoissonMod, only: Padepoisson 
    use sgsmod, only: sgs
    use wallmodelMod, only: wallmodel
    use numerics
    use cd06staggstuff, only: cd06stagg
    use cf90stuff, only: cf90
    use PadePoissonMod, only: Padepoisson 

    implicit none

    private
    public :: igridWallM 

    complex(rkind), parameter :: zeroC = zero + imi*zero 
    integer, parameter :: no_slip = 1, slip = 2



    ! Allow non-zero value (isEven) 
    logical :: topBC_u = .true.  , topBC_v = .true. 
    logical :: botBC_u = .false. , botBC_v = .false. 
    logical, parameter :: topBC_w = .false. , botBC_w = .false. 
    integer :: ierr 

    type, extends(grid) :: igridWallM
        
        character(clen) :: inputDir

        type(decomp_info), allocatable :: gpC, gpE
        type(decomp_info), pointer :: Sp_gpC, Sp_gpE
        type(spectral), allocatable :: spectE, spectC
        type(staggOps), allocatable :: Ops, OpsPP
        type(sgs), allocatable :: sgsmodel
        type(wallmodel), allocatable :: moengWall

        real(rkind), dimension(:,:,:,:), allocatable :: PfieldsC
        real(rkind), dimension(:,:,:,:), allocatable :: PfieldsE
        type(cd06stagg), allocatable :: derW, derWW, derOO, derEE
        type(cf90),      allocatable :: filzE, filzC

        complex(rkind), dimension(:,:,:,:), allocatable :: SfieldsC
        complex(rkind), dimension(:,:,:,:), allocatable :: SfieldsE


        type(poisson), allocatable :: poiss
        type(padepoisson), allocatable :: padepoiss
        real(rkind), dimension(:,:,:), allocatable :: divergence

        real(rkind), dimension(:,:,:), pointer :: u, v, wC, w, uE, vE
        real(rkind), dimension(:,:,:), pointer :: ox,oy,oz
        complex(rkind), dimension(:,:,:), pointer :: uhat, vhat, whatC, what
        real(rkind), dimension(:,:,:), pointer :: T, TE
        complex(rkind), dimension(:,:,:), pointer :: That, TEhat, T_rhs, T_Orhs

        complex(rkind), dimension(:,:,:), allocatable :: uBase, Tbase, dTdxH, dTdyH, dTdzH
        real(rkind), dimension(:,:,:), allocatable :: dTdxC, dTdyC, dTdzE

        real(rkind), dimension(:,:,:,:), allocatable, public :: rbuffxC, rbuffyC, rbuffzC
        real(rkind), dimension(:,:,:,:), allocatable :: rbuffxE, rbuffyE, rbuffzE
        
        complex(rkind), dimension(:,:,:,:), allocatable :: cbuffyC, cbuffzC
        complex(rkind), dimension(:,:,:,:), allocatable :: cbuffyE, cbuffzE

        complex(rkind), dimension(:,:,:,:), allocatable :: rhsC, rhsE, OrhsC, OrhsE 
        real(rkind), dimension(:,:,:,:), allocatable :: duidxjC, duidxjE 
        complex(rkind), dimension(:,:,:,:), allocatable :: duidxjChat
        complex(rkind), dimension(:,:,:), pointer:: u_rhs, v_rhs, wC_rhs, w_rhs 
        complex(rkind), dimension(:,:,:), pointer:: u_Orhs, v_Orhs, w_Orhs

        real(rkind), dimension(:,:,:), allocatable :: rDampC, rDampE         
        real(rkind) :: Re, Gx, Gy, Gz, dtby2, meanfact, Tref
        complex(rkind), dimension(:,:,:), allocatable :: GxHat 
        real(rkind) :: Ro = 1.d5, gravity_nd = 9.8d0, Fr = 1000.d0

        integer :: nxZ, nyZ
        
        integer :: runID, t_start_planeDump, t_stop_planeDump, t_planeDump, t_DivergenceCheck
        logical :: useCoriolis = .true. , isStratified = .false., useSponge = .false. 
        logical :: useExtraForcing = .false., useGeostrophicForcing = .false. 
        logical :: useSGS = .false. 
        logical :: UseDealiasFilterVert = .false.
        logical :: useDynamicProcedure 
        logical :: useCFL = .false.  
        logical :: dumpPlanes = .false.

        complex(rkind), dimension(:,:,:), allocatable :: dPf_dxhat

        real(rkind) :: max_nuSGS, invObLength, Tsurf, dTsurf_dt

        real(rkind) :: z0, ustar = zero, Umn, Vmn, Uspmn, dtOld, dtRat, Tmn, wTh_surf
        real(rkind), dimension(:,:,:), allocatable :: filteredSpeedSq
        integer :: wallMType, botBC_Temp 

        ! Statistics to compute 
        real(rkind), dimension(:,:), allocatable :: zStats2dump, runningSum, TemporalMnNOW
        real(rkind), dimension(:), pointer :: u_mean, v_mean, w_mean, uu_mean, uv_mean, uw_mean, vv_mean, vw_mean, ww_mean
        real(rkind), dimension(:), pointer :: tau11_mean, tau12_mean, tau13_mean, tau22_mean, tau23_mean, tau33_mean
        real(rkind), dimension(:), pointer :: S11_mean, S12_mean, S13_mean, S22_mean, S23_mean, S33_mean
        real(rkind), dimension(:), pointer :: viscdissp, sgsdissp, sgscoeff_mean, PhiM
        integer :: tidSUM, tid_StatsDump, tid_compStats,tSimStartStats
        logical :: normByustar

        ! Pointers linked to SGS stuff
        real(rkind), dimension(:,:,:,:), pointer :: tauSGS_ij
        real(rkind), dimension(:,:,:)  , pointer :: nu_SGS, tau13, tau23
        real(rkind), dimension(:,:,:)  , pointer :: c_SGS 
        
        integer, dimension(:), allocatable :: xplanes, yplanes, zplanes
        ! Note that c_SGS is linked to a variable that is constant along & 
        ! i, j but is still stored as a full 3 rank array. This is mostly done to
        ! make it convenient us to later do transposes or to compute Sij.


        contains
            procedure :: init
            procedure :: init_stats
            procedure :: destroy
            procedure :: printDivergence 
            procedure :: AdamsBashforth
            procedure :: getMaxKE
            procedure, private :: interp_primitiveVars
            procedure, private :: compute_duidxj
            procedure, private :: compute_dTdxi
            procedure, private :: addNonLinearTerm_Rot
            procedure, private :: AddBuoyancyTerm
            procedure, private :: addCoriolisTerm
            procedure, private :: addSponge
            procedure, private :: addExtraForcingTerm 
            procedure, private :: compute_and_bcast_surface_Mn
            procedure, private :: dumpRestartFile
            procedure, private :: readRestartFile
            procedure, private :: compute_z_mean 
            procedure, private :: compute_z_fluct
            procedure, private :: compute_deltaT
            procedure, private :: getfilteredSpeedSqAtWall
            procedure, private :: dump_stats
            procedure, private :: compute_stats 
            procedure, private :: getSurfaceQuantities 
            procedure, private :: ApplyCompactFilter 
            procedure, private :: addNonLinearTerm_skewSymm
            procedure          :: finalize_stats
            procedure, private :: dump_planes
            procedure          :: dumpFullField 
    end type

contains 

    subroutine init(this,inputfile)
        class(igridWallM), intent(inout), target :: this        
        character(len=clen), intent(in) :: inputfile 
        character(len=clen) :: outputdir, inputdir
        integer :: nx, ny, nz, prow = 0, pcol = 0, ioUnit, nsteps = -1, topWall = slip, SGSModelID = 1
        integer :: tid_StatsDump =10000, tid_compStats = 10000,  WallMType = 0, t_planeDump = 1000
        integer :: runID = 0,  t_dataDump = 99999, t_restartDump = 99999, t_stop_planeDump = 1
        integer :: restartFile_TID = 1, ioType = 0, restartFile_RID =1, t_start_planeDump = 1, botBC_Temp = 0
        real(rkind) :: dt=-one,tstop=one,CFL =-one,tSimStartStats=100.d0,dpfdy=zero,dPfdz=zero,ztop,ncWall=1.d0
        real(rkind) :: Pr = 0.7_rkind, Re = 8000._rkind, Ro = 1000._rkind,dpFdx = zero, z0 = 1.d-4, Cs = 0.17d0
        real(rkind) :: SpongeTscale = 50._rkind, zstSponge = 0.8_rkind, Fr = 1000.d0,Gx=0.d0,Gy=0.d0,Gz=0.d0
        logical :: useRestartFile=.false.,isInviscid=.false.,useCoriolis = .true. 
        logical :: isStratified=.false.,dumpPlanes = .false., useSGSclipping = .true.,useExtraForcing = .false.
        logical :: useSGS = .true., useDynamicProcedure = .false., useSpongeLayer=.false.
        logical :: useGeostrophicForcing = .false., useVerticalTfilter = .false., useWallDamping = .true. 
        real(rkind), dimension(:,:,:), pointer :: zinZ, zinY, zEinY, zEinZ
        integer :: AdvectionTerm = 1, NumericalSchemeVert = 0, t_DivergenceCheck = 10
        logical :: normStatsByUstar=.false., ComputeStokesPressure = .false., UseDealiasFilterVert = .false.

        namelist /INPUT/ nx, ny, nz, tstop, dt, CFL, nsteps, inputdir, outputdir, prow, pcol, &
                         useRestartFile, restartFile_TID, restartFile_RID 
        namelist /IO/ t_restartDump, t_dataDump, ioType, dumpPlanes, runID, &
                        t_planeDump, t_stop_planeDump, t_start_planeDump 
        namelist /STATS/ tid_StatsDump, tid_compStats, tSimStartStats, normStatsByUstar
        namelist /PHYSICS/isInviscid,useCoriolis,useExtraForcing,isStratified,Re,Ro,Pr,Fr, &
                          useGeostrophicForcing, Gx, Gy, Gz, dpFdx, dpFdy, dpFdz
        namelist /BCs/ topWall, useSpongeLayer, zstSponge, SpongeTScale, botBC_Temp
        namelist /LES/ useSGS, useDynamicProcedure, useSGSclipping, SGSmodelID, useVerticalTfilter, &
                        useWallDamping, ncWall, Cs 
        namelist /WALLMODEL/ z0, wallMType 
        namelist /NUMERICS/ AdvectionTerm, ComputeStokesPressure, NumericalSchemeVert, &
                            UseDealiasFilterVert, t_DivergenceCheck


        ! STEP 1: READ INPUT 
        ioUnit = 11
        open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
        read(unit=ioUnit, NML=INPUT)
        read(unit=ioUnit, NML=NUMERICS)
        read(unit=ioUnit, NML=IO)
        read(unit=ioUnit, NML=STATS)
        read(unit=ioUnit, NML=PHYSICS)
        read(unit=ioUnit, NML=BCs)
        read(unit=ioUnit, NML=LES)
        read(unit=ioUnit, NML=WALLMODEL)
        close(ioUnit)
      
        this%nx = nx; this%ny = ny; this%nz = nz; this%meanfact = one/(real(nx,rkind)*real(ny,rkind)); 
        this%dt = dt; this%dtby2 = dt/two ; this%z0 = z0 ; this%Re = Re; this%useSponge = useSpongeLayer
        this%outputdir = outputdir; this%inputdir = inputdir; this%isStratified = isStratified 
        this%WallMtype = WallMType; this%runID = runID; this%tstop = tstop; this%t_dataDump = t_dataDump
        this%CFL = CFL; this%dumpPlanes = dumpPlanes; this%useGeostrophicForcing = useGeostrophicForcing
        if (this%CFL > zero) this%useCFL = .true. 
        if ((this%CFL < zero) .and. (this%dt < zero)) then
            call GracefulExit("Both CFL and dt cannot be negative. Have you &
            & specified either one of these in the input file?", 124)
        end if 
        this%t_restartDump = t_restartDump; this%tid_statsDump = tid_statsDump
        this%useCoriolis = useCoriolis; this%tSimStartStats = tSimStartStats
        this%tid_compStats = tid_compStats; this%useExtraForcing = useExtraForcing; this%useSGS = useSGS 
        this%useDynamicProcedure = useDynamicProcedure; this%UseDealiasFilterVert = UseDealiasFilterVert
        this%Gx = Gx; this%Gy = Gy; this%Gz = Gz; this%Fr = Fr; this%gravity_nd = one/(this%Fr**2)
        this%t_start_planeDump = t_start_planeDump; this%t_stop_planeDump = t_stop_planeDump
        this%t_planeDump = t_planeDump; this%BotBC_temp = BotBC_temp; this%Ro = Ro
        this%normByustar = normStatsByUstar; this%t_DivergenceCheck = t_DivergenceCheck
        


        ! STEP 2: ALLOCATE DECOMPOSITIONS
        allocate(this%gpC); allocate(this%gpE)
        call decomp_2d_init(nx, ny, nz, prow, pcol)
        call get_decomp_info(this%gpC)
        call decomp_info_init(nx,ny,nz+1,this%gpE)
       
        ! Set Top and Bottom BCs
        if (topWall == slip) then
            topBC_u = .true.; topBC_v = .true.
            call message(1, "TopWall BC set to: SLIP")
        elseif (topWall == no_slip) then
            topBC_u = .false.; topBC_v = .false.
            call message(1, "TopWall BC set to: NO_SLIP")
        elseif (topWall == 3) then
            call GracefulExit("Wall model for the top wall is not currently supported", 321)
        else
            call message("WARNING: No Top BCs provided. Using defaults found in igridWallM.F90")
        end if 



        ! Set numerics
        select case(NumericalSchemeVert)
        case(0)
            useCompactFD = .false.
        case(1)
            useCompactFD = .true.
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
        
        ! STEP 4: ALLOCATE/INITIALIZE THE SPECTRAL DERIVED TYPES
        allocate(this%spectC)
        call this%spectC%init("x", nx, ny, nz, this%dx, this%dy, this%dz, &
                "four", this%filter_x, 2 , .false.)
        allocate(this%spectE)
        call this%spectE%init("x", nx, ny, nz+1, this%dx, this%dy, this%dz, &
                "four", this%filter_x, 2 , .false.)
        this%sp_gpC => this%spectC%spectdecomp
        this%sp_gpE => this%spectE%spectdecomp


        ! STEP 5: ALLOCATE/INITIALIZE THE OPERATORS DERIVED TYPE
        if (useCompactFD) then
            allocate(this%derEE, this%derOO, this%derW, this%derWW) 
            call this%derEE%init( this%gpC%zsz(3), this%dz, isTopEven = .true., isBotEven = .true., & 
                             isTopSided = .false., isBotSided = .true.) 
            call this%derOO%init( this%gpC%zsz(3), this%dz, isTopEven = .false., isBotEven = .false., & 
                             isTopSided = .false., isBotSided = .true.) 
            call this%derW%init( this%gpC%zsz(3),  this%dz, isTopEven = .false., isBotEven = .false., & 
                             isTopSided = .false., isBotSided = .false.) 
            call this%derWW%init( this%gpC%zsz(3),  this%dz, isTopEven = .true., isBotEven = .true., & 
                             isTopSided = .false., isBotSided = .false.) 
        else
            allocate(this%Ops)
            call this%Ops%init(this%gpC,this%gpE,0,this%dx,this%dy,this%dz,this%spectC%spectdecomp, &
                    this%spectE%spectdecomp, .false., .false.)
        end if        
        allocate(this%OpsPP)
        call this%OpsPP%init(this%gpC,this%gpE,0,this%dx,this%dy,this%dz,this%spectC%spectdecomp, &
                    this%spectE%spectdecomp, .false., .false.)
        
        if (this%UseDealiasFilterVert) then
            allocate(this%filzC, this%filzE)
            ierr = this%filzC%init(nz  , .false.)
            ierr = this%filzE%init(nz+1, .false.)
        end if


 
        ! STEP 6: ALLOCATE MEMORY FOR FIELD ARRAYS
        if (this%isStratified) then
            allocate(this%PfieldsC(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3),7))
            allocate(this%PfieldsE(this%gpE%xsz(1),this%gpE%xsz(2),this%gpE%xsz(3),4))
            allocate(this%dTdzE(this%gpE%xsz(1),this%gpE%xsz(2),this%gpE%xsz(3)))
            allocate(this%dTdxC(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)))
            allocate(this%dTdyC(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)))
            call this%spectC%alloc_r2c_out(this%SfieldsC,4)
            call this%spectC%alloc_r2c_out(this%dTdxH)
            call this%spectC%alloc_r2c_out(this%dTdyH)
            call this%spectE%alloc_r2c_out(this%dTdzH)
            call this%spectC%alloc_r2c_out(this%rhsC,3); call this%spectC%alloc_r2c_out(this%OrhsC,3)
            call this%spectE%alloc_r2c_out(this%SfieldsE,2)
        else
            allocate(this%PfieldsC(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3),6))
            allocate(this%PfieldsE(this%gpE%xsz(1),this%gpE%xsz(2),this%gpE%xsz(3),3))
            call this%spectC%alloc_r2c_out(this%SfieldsC,3)
            call this%spectC%alloc_r2c_out(this%rhsC,2); call this%spectC%alloc_r2c_out(this%OrhsC,2)
            call this%spectE%alloc_r2c_out(this%SfieldsE,1)
        end if 
        allocate(this%divergence(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)))
        allocate(this%duidxjC(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3),9))
        allocate(this%duidxjE(this%gpE%xsz(1),this%gpE%xsz(2),this%gpE%xsz(3),4))
        call this%spectC%alloc_r2c_out(this%duidxjChat,9)
        call this%spectE%alloc_r2c_out(this%rhsE,1); call this%spectE%alloc_r2c_out(this%OrhsE,1)
        
        this%u => this%PfieldsC(:,:,:,1) ; this%v => this%PfieldsC(:,:,:,2) ; this%wC => this%PfieldsC(:,:,:,3) 
        this%w => this%PfieldsE(:,:,:,1) ; this%uE => this%PfieldsE(:,:,:,2) ; this%vE => this%PfieldsE(:,:,:,3) 
        
        this%uhat => this%SfieldsC(:,:,:,1); this%vhat => this%SfieldsC(:,:,:,2); 
        this%whatC => this%SfieldsC(:,:,:,3); this%what => this%SfieldsE(:,:,:,1)

        this%ox => this%PfieldsC(:,:,:,4); this%oy => this%PfieldsC(:,:,:,5); this%oz => this%PfieldsC(:,:,:,6)

        this%u_rhs => this%rhsC(:,:,:,1); this%v_rhs => this%rhsC(:,:,:,2); this%w_rhs => this%rhsE(:,:,:,1)

        this%u_Orhs => this%OrhsC(:,:,:,1); this%v_Orhs => this%OrhsC(:,:,:,2); this%w_Orhs => this%OrhsE(:,:,:,1)

        if (this%isStratified) then
            this%T => this%PfieldsC(:,:,:,7); this%That => this%SfieldsC(:,:,:,4)
            this%TE => this%PfieldsE(:,:,:,4); this%T_rhs => this%rhsC(:,:,:,3)
            this%T_Orhs => this%OrhsC(:,:,:,3); this%TEhat => this%SfieldsE(:,:,:,2)
        end if

        allocate(this%cbuffyC(this%sp_gpC%ysz(1),this%sp_gpC%ysz(2),this%sp_gpC%ysz(3),2))
        allocate(this%cbuffyE(this%sp_gpE%ysz(1),this%sp_gpE%ysz(2),this%sp_gpE%ysz(3),2))
        
        allocate(this%cbuffzC(this%sp_gpC%zsz(1),this%sp_gpC%zsz(2),this%sp_gpC%zsz(3),2))
        allocate(this%cbuffzE(this%sp_gpE%zsz(1),this%sp_gpE%zsz(2),this%sp_gpE%zsz(3),2))

        allocate(this%rbuffxC(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3),2))
        allocate(this%rbuffyC(this%gpC%ysz(1),this%gpC%ysz(2),this%gpC%ysz(3),2))
        allocate(this%rbuffzC(this%gpC%zsz(1),this%gpC%zsz(2),this%gpC%zsz(3),4))

        allocate(this%rbuffxE(this%gpE%xsz(1),this%gpE%xsz(2),this%gpE%xsz(3),2))
        allocate(this%rbuffyE(this%gpE%ysz(1),this%gpE%ysz(2),this%gpE%ysz(3),2))
        allocate(this%rbuffzE(this%gpE%zsz(1),this%gpE%zsz(2),this%gpE%zsz(3),4))
        allocate(this%filteredSpeedSq(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)))
        this%nxZ = size(this%cbuffzE,1); this%nyZ = size(this%cbuffzE,2)

        ! STEP 6: ALLOCATE/INITIALIZE THE POISSON DERIVED TYPE
        if (useCompactFD) then
            allocate(this%padepoiss)
            call this%padepoiss%init(this%dx, this%dy, this%dz, this%spectC, this%spectE, this%derW) 
        else    
            allocate(this%poiss)
            call this%poiss%init(this%spectC,.false.,this%dx,this%dy,this%dz,this%Ops,this%spectE, .true.)  
        end if 
               
        ! STEP 7: INITIALIZE THE FIELDS
        if (useRestartFile) then
            call this%readRestartFile(restartfile_TID, restartfile_RID)
            this%step = restartfile_TID
        else 
            call initfields_wallM(this%gpC, this%gpE, inputfile, this%mesh, this%PfieldsC, this%PfieldsE)! <-- this procedure is part of user defined HOOKS
            this%step = 0
            this%tsim = zero
            call this%dumpRestartfile()
        end if 
       
        if (botBC_Temp == 0) then
            call setDirichletBC_Temp(inputfile, this%Tsurf, this%dTsurf_dt)
        else
            call GraceFulExit("Only Dirichlet BC supported for Temperature at &
                & this time. Set botBC_Temp = 0",341)        
        end if 

        call this%spectC%fft(this%u,this%uhat)   
        call this%spectC%fft(this%v,this%vhat)   
        call this%spectE%fft(this%w,this%what)   
        if (this%isStratified) call this%spectC%fft(this%T,this%That)   

        ! Dealias and filter before projection
        call this%spectC%dealias(this%uhat)
        call this%spectC%dealias(this%vhat)
        call this%spectE%dealias(this%what)
        if (this%isStratified) call this%spectC%dealias(this%That)
        if (this%UseDealiasFilterVert) then
            call this%ApplyCompactFilter()
        end if


        ! Pressure projection
        if (useCompactFD) then
            call this%padepoiss%PressureProjection(this%uhat,this%vhat,this%what)
            call this%padepoiss%DivergenceCheck(this%uhat, this%vhat, this%what, this%divergence)
        else
            call this%poiss%PressureProjNP(this%uhat,this%vhat,this%what)
            call this%poiss%DivergenceCheck(this%uhat, this%vhat, this%what, this%divergence)
        end if 

        ! Take it back to physical fields
        call this%spectC%ifft(this%uhat,this%u)
        call this%spectC%ifft(this%vhat,this%v)
        call this%spectE%ifft(this%what,this%w)
        if (this%isStratified) call this%spectC%ifft(this%That,this%T)

        ! STEP 8: Interpolate the cell center values of w
        call this%compute_and_bcast_surface_Mn()
        call this%interp_PrimitiveVars()
        call message(1,"Max KE:",P_MAXVAL(this%getMaxKE()))
     
        ! STEP 9: Compute duidxj
        call this%compute_duidxj()
        if (this%isStratified) call this%compute_dTdxi() 

        ! STEP 10a: Compute Coriolis Term
        if (this%useCoriolis) then
            call message(0, "Turning on Coriolis with Geostrophic Forcing")
            call message(1, "Geostrophic Velocity:", this%Gx) 
            call message(1, "Rossby Number       :", this%Ro) 
            call this%spectC%alloc_r2c_out(this%GxHat)
            this%rbuffxC(:,:,:,1) = this%Gx
            call this%spectC%fft(this%rbuffxC(:,:,:,1),this%Gxhat)
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
            allocate(this%SGSmodel)
            call this%SGSmodel%init(SGSModelID, this%spectC, this%spectE, this%gpC, this%gpE, this%dx, & 
                this%dy, this%dz, useDynamicProcedure, useSGSclipping, this%mesh(:,:,:,3), this%z0, &
                .true., WallMType, useVerticalTfilter, Pr, useWallDamping, nCWall, Cs, ComputeStokesPressure )
            call this%sgsModel%link_pointers(this%nu_SGS, this%c_SGS, this%tauSGS_ij, this%tau13, this%tau23)
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
            ztop = zEinZ(1,1,this%nz+1)
            call transpose_z_to_y(zEinZ,zEinY,this%gpE)
            this%RdampC = (one/SpongeTscale) * (one - cos(pi*(zinY - zstSponge) /(zTop - zstSponge)))/two
            this%RdampE = (one/SpongeTscale) * (one - cos(pi*(zEinY - zstSponge)/(zTop - zstSponge)))/two
            where (zEinY < zstSponge) 
                this%RdampE = zero
            end where
            where (zinY < zstSponge) 
                this%RdampC = zero
            end where
            call message(0,"Sponge Layer initialized successfully")
        end if 


        ! STEP 12: Set visualization planes for io
        call set_planes_io(this%xplanes, this%yplanes, this%zplanes)


        ! STEP 13: Compute the timestep
        call this%compute_deltaT()
        this%dtOld = this%dt
        this%dtRat = one 

        ! Final Step: Safeguard against unfinished procedures
        if ((useCompactFD).and.(.not.useSkewSymm)) then
            call GracefulExit("You must solve in skew symmetric form if you use CD06",54)
        end if 
        if ((useCompactFD).and.(ComputeStokesPressure)) then
            call GracefulExit("You cannot compute Stokes Pressure if you use CD06",54)
        end if 
        if ((useCompactFD).and.(this%isStratified)) then
            call GracefulExit("Stratified case is not currently supported using CD06",54)
        end if 

        call message("IGRID initialized successfully!")
        call message("===========================================================")


    end subroutine


    subroutine compute_deltaT(this)
        use reductions, only: p_minval
        class(igridWallM), intent(inout), target :: this
        real(rkind) :: t1, t2, t3, tmin


        if (this%useCFL) then
            t1 = p_maxval(maxval(this%u))
            t1 = this%dx/(t1 + 1d-13)

            t2 = p_maxval(maxval(this%v))
            t2 = this%dy/(t2 + 1d-13)
            
            t3 = p_maxval(maxval(this%w))
            t3 = this%dz/(t3 + 1d-13)
            
            tmin = min(t1,t2,t3)
            this%dt = this%CFL*tmin
        end if 


    end subroutine


    function getMaxKE(this) result(maxKE)
        class(igridWallM), intent(inout) :: this
        real(rkind)  :: maxKE

        this%rbuffxC(:,:,:,1) = this%u**2 + this%v**2 + this%wC**2
        maxKE = half*p_maxval(maxval(this%rbuffxC))

    end function

    subroutine interp_PrimitiveVars(this)
        class(igridWallM), intent(inout), target :: this
        complex(rkind), dimension(:,:,:), pointer :: ybuffC, ybuffE, zbuffC, zbuffE

        ybuffE => this%cbuffyE(:,:,:,1)
        zbuffE => this%cbuffzE(:,:,:,1)
        zbuffC => this%cbuffzC(:,:,:,1)
        ybuffC => this%cbuffyC(:,:,:,1)


        ! Step 1: Interpolate w -> wC
        call transpose_y_to_z(this%what,zbuffE,this%sp_gpE)
        if (useCompactFD) then
            call this%derW%InterpZ_E2C(zbuffE,zbuffC,size(zbuffE,1),size(zbuffE,2))
        else
            call this%Ops%InterpZ_Edge2Cell(zbuffE,zbuffC)
        end if
        call transpose_z_to_y(zbuffC,this%whatC,this%sp_gpC)
        call this%spectC%ifft(this%whatC,this%wC)

        ! Step 2: Interpolate u -> uE
        call transpose_y_to_z(this%uhat,zbuffC,this%sp_gpC)
        if (useCompactFD) then
            call this%derEE%InterpZ_C2E(zbuffC,zbuffE,size(zbuffC,1),size(zbuffC,2))
        else
            call this%Ops%InterpZ_Cell2Edge(zbuffC,zbuffE,zeroC,zeroC)
            zbuffE(:,:,this%nz + 1) = zbuffC(:,:,this%nz)
        end if 
        zbuffE(:,:,1) = (three/two)*zbuffC(:,:,1) - half*zbuffC(:,:,2)

        call transpose_z_to_y(zbuffE,ybuffE,this%sp_gpE)
        call this%spectE%ifft(ybuffE,this%uE)
        
        ! Step 3: Interpolate v -> vE
        call transpose_y_to_z(this%vhat,zbuffC,this%sp_gpC)
        if (useCompactFD) then
            call this%derEE%InterpZ_C2E(zbuffC,zbuffE,size(zbuffC,1),size(zbuffC,2))
        else
            call this%Ops%InterpZ_Cell2Edge(zbuffC,zbuffE,zeroC,zeroC)
            zbuffE(:,:,this%nz + 1) = zbuffC(:,:,this%nz)
        end if 
        
        zbuffE(:,:,1) = (three/two)*zbuffC(:,:,1) - half*zbuffC(:,:,2)
        
        call transpose_z_to_y(zbuffE,ybuffE,this%sp_gpE)
        call this%spectE%ifft(ybuffE,this%vE)

        ! Step 4: Interpolate T
        if (this%isStratified) then
            call transpose_y_to_z(this%That,zbuffC,this%sp_gpC)
            call this%Ops%InterpZ_Cell2Edge(zbuffC,zbuffE,zeroC,zeroC)
            zbuffE(:,:,1) = (three/two)*zbuffC(:,:,1) - half*zbuffC(:,:,2)
            zbuffE(:,:,this%nz + 1) = two*zbuffC(:,:,this%nz) - zbuffE(:,:,this%nz)
            call transpose_z_to_y(zbuffE,this%TEhat,this%sp_gpE)
            call this%spectE%ifft(this%TEhat,this%TE)
        end if 
    end subroutine


    subroutine printDivergence(this)
        class(igridWallM), intent(inout) :: this
        if (useCompactFD) then
            call this%padepoiss%DivergenceCheck(this%uhat, this%vhat, this%what, this%divergence)
        else
            call this%poiss%DivergenceCheck(this%uhat, this%vhat, this%what, this%divergence)
        end if 
    end subroutine 


    subroutine destroy(this)
        class(igridWallM), intent(inout) :: this
        
        nullify(this%u, this%uhat, this%v, this%vhat, this%w, this%what, this%wC)
        deallocate(this%PfieldsC, this%PfieldsE, this%SfieldsC, this%SfieldsE)
        nullify(this%u_rhs, this%v_rhs, this%w_rhs)
        deallocate(this%rhsC, this%rhsE, this%OrhsC, this%OrhsE)
        deallocate(this%duidxjC, this%duidxjChat)
        call this%spectC%destroy()
        call this%spectE%destroy()
        deallocate(this%spectC, this%spectE)
        nullify(this%nu_SGS, this%c_SGS, this%tauSGS_ij)
        call this%sgsModel%destroy()
        deallocate(this%sgsModel)
    end subroutine

    subroutine addNonLinearTerm_Rot(this)
        class(igridWallM), intent(inout), target :: this
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

        dwdx => this%duidxjE(:,:,:,1); dwdy => this%duidxjE(:,:,:,2);
        dudz => this%duidxjE(:,:,:,3); dvdz => this%duidxjE(:,:,:,4);

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
        call this%Ops%InterpZ_Edge2Cell(tzE,tzC)
        call transpose_z_to_y(tzC,this%u_rhs, this%sp_gpC)
        this%u_rhs = this%u_rhs + fT1C


        T1C = dudy - dvdx
        T1C = T1C*this%u
        call this%spectC%fft(T1C,fT1C)
        T2E = dwdy - dvdz
        T2E = T2E*this%w
        call this%spectE%fft(T2E,fT2E)
        call transpose_y_to_z(fT2E,tzE, this%sp_gpE)
        call this%Ops%InterpZ_Edge2Cell(tzE,tzC)
        call transpose_z_to_y(tzC,this%v_rhs, this%sp_gpC)
        this%v_rhs = this%v_rhs + fT1C

        T1E = dudz - dwdx
        T1E = T1E*this%uE
        T2E = dvdz - dwdy
        T2E = T2E*this%vE
        T1E = T1E + T2E
        call this%spectE%fft(T1E,this%w_rhs)

        if (this%isStratified) then
            T1C = -this%u*this%dTdxC 
            T2C = -this%v*this%dTdyC
            T1C = T1C + T2C
            call this%spectC%fft(T1c,fT1C)
            T1E = -this%w*this%dTdzE
            call this%spectE%fft(T1E,fT1E)
            call transpose_y_to_z(fT2E,tzE, this%sp_gpE)
            call this%Ops%InterpZ_Edge2Cell(tzE,tzC)
            call transpose_z_to_y(tzC,this%T_rhs, this%sp_gpC)
            this%T_rhs = this%T_rhs + fT1C
        end if

    end subroutine

    subroutine addNonLinearTerm_skewSymm(this)
        class(igridWallM), intent(inout), target :: this
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

        dwdx => this%duidxjE(:,:,:,1); dwdy => this%duidxjE(:,:,:,2);
        dudz => this%duidxjE(:,:,:,3); dvdz => this%duidxjE(:,:,:,4);

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
        if (useCompactFD) then
            call this%derW%InterpZ_E2C(tzE,tzC,size(tzE,1),size(tzE,2))
        else
            call this%Ops%InterpZ_Edge2Cell(tzE,tzC)
        end if 
        call transpose_z_to_y(tzC,this%u_rhs, this%sp_gpC)
        this%u_rhs = this%u_rhs + fT1C
        
        T1C = dvdx*this%u
        T2C = dvdy*this%v
        T1C = T1C + T2C
        T1E = dvdz*this%w
        call this%spectC%fft(T1C,fT1C)
        call this%spectE%fft(T1E,fT1E)
        call transpose_y_to_z(fT1E,tzE, this%sp_gpE)
        if (useCompactFD) then
            call this%derW%InterpZ_E2C(tzE,tzC,size(tzE,1),size(tzE,2))
        else
            call this%Ops%InterpZ_Edge2Cell(tzE,tzC)
        end if 
        call transpose_z_to_y(tzC,this%v_rhs, this%sp_gpC)
        this%v_rhs = this%v_rhs + fT1C
        
        T1E = dwdx*this%uE
        T2E = dwdy*this%vE
        T2E = T1E + T2E
        call this%spectE%fft(T2E,fT2E)
        T1C = dwdz*this%wC
        call this%spectC%fft(T1C,fT1C)
        call transpose_y_to_z(fT1C,tzC, this%sp_gpC)
        if (useCompactFD) then
            call this%derW%InterpZ_C2E(tzC,tzE,size(tzC,1),size(tzC,2))
        else
            call this%Ops%InterpZ_Cell2Edge(tzC,tzE,zeroC,zeroC)
        end if 
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
        if (useCompactFD) then
            call this%derWW%ddz_C2E(tzC,tzE,size(tzC,1),size(tzC,2))
        else
            call this%Ops%ddz_C2E(tzC,tzE,.true.,.true.)
        end if
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
        if (useCompactFD) then
            call this%derW%ddz_E2C(tzE,tzC,size(tzE,1),size(tzE,2))
        else
            call this%Ops%ddz_E2C(tzE,tzC)
        end if 
        call transpose_z_to_y(tzC,fT1C,this%sp_gpC)
        this%u_rhs = this%u_rhs + fT1C
        
        call this%spectE%mtimes_ik1_ip(fT1E)
        this%w_rhs = this%w_rhs + fT1E


        T1E = this%vE*this%w
        call this%spectE%fft(T1E,fT1E)
        call transpose_y_to_z(fT1E,TzE,this%sp_gpE)
        if (useCompactFD) then
            call this%derW%ddz_E2C(tzE,tzC,size(tzE,1),size(tzE,2))
        else
            call this%Ops%ddz_E2C(tzE,tzC)
        end if 
        call transpose_z_to_y(tzC,fT1C,this%sp_gpC)
        this%v_rhs = this%v_rhs + fT1C

        call this%spectE%mtimes_ik2_ip(fT1E)
        this%w_rhs = this%w_rhs + fT1E

        this%u_rhs = -half*this%u_rhs
        this%v_rhs = -half*this%v_rhs
        this%w_rhs = -half*this%w_rhs

    end subroutine

    subroutine addCoriolisTerm(this)
        class(igridWallM), intent(inout) :: this
        ! u equation 
        this%u_rhs = this%u_rhs + this%vhat/this%Ro
        ! v equation 
        this%v_rhs = this%v_rhs +  (this%GxHat - this%uhat)/this%Ro
    end subroutine  

    subroutine addSponge(this)
        class(igridWallM), intent(inout), target :: this
        complex(rkind), dimension(:,:,:), pointer :: deviationC
        deviationC => this%cbuffyC(:,:,:,1)
        
        deviationC = this%uhat - this%ubase
        this%u_rhs = this%u_rhs - (this%RdampC/this%dt)*deviationC

        this%v_rhs = this%v_rhs - (this%RdampC/this%dt)*this%vhat ! base value for v is zero
        
        this%w_rhs = this%w_rhs - (this%RdampE/this%dt)*this%what ! base value for w is zero  

        deviationC = this%That - this%Tbase
        this%T_rhs = this%T_rhs - (this%RdampC/this%dt)*deviationC

    end subroutine

    subroutine addExtraForcingTerm(this)
        class(igridWallM), intent(inout) :: this
        !if (this%spectC%carryingZeroK) then
        !    this%dpF_dxhat(1,1,:) = cmplx(this%ustar*this%ustar*this%nx*this%ny,zero)
        !end if
        this%u_rhs = this%u_rhs + this%dpF_dxhat
    end subroutine

    subroutine AddBuoyancyTerm(this)
        class(igridWallM), intent(inout), target :: this
        complex(rkind), dimension(:,:,:), pointer :: fT1E 
   
        fT1E => this%cbuffyE(:,:,:,1)
        fT1E = this%TEhat/(this%Fr*this%Fr)
        if (this%spectE%carryingZeroK) then
            fT1E(1,1,:) = zero
        end if 
        this%w_rhs = this%w_rhs + fT1E 
    end subroutine

    subroutine AdamsBashforth(this)
        class(igridWallM), intent(inout) :: this
        real(rkind) :: abf1, abf2

        ! Step 0: Compute TimeStep 
        call this%compute_deltaT
        this%dtRat = this%dt/this%dtOld

        ! Step 1: Non Linear Term 
        if (useSkewSymm) then
            call this%addNonLinearTerm_skewSymm()
        else
            call this%AddNonLinearTerm_Rot()
        end if 

        ! Step 2: Coriolis Term
        if (this%useCoriolis) then
            call this%AddCoriolisTerm()
        end if 
      
        ! Step 3: Extra Forcing 
        if (this%useExtraForcing) then
            call this%addExtraForcingTerm()
        end if 

        ! Step 4: Buoyancy Term
        if (this%isStratified) then
            call this%addBuoyancyTerm()
        end if 

        
        ! Step 4: SGS Viscous Term
        if (this%useSGS) then
            call this%SGSmodel%getRHS_SGS_WallM(this%duidxjC, this%duidxjE        , this%duidxjChat ,& 
                                                this%u_rhs  , this%v_rhs          , this%w_rhs      ,&
                                                this%uhat   , this%vhat           , this%whatC      ,&
                                                this%u      , this%v              , this%wC         ,&
                                                this%ustar  , this%Umn            , this%Vmn        ,&
                                                this%Uspmn  , this%filteredSpeedSq, this%InvObLength,&
                                                this%max_nuSGS)

            !! IMPORTANT: duidxjC, u, v and wC are all corrupted if SGS was initialized to use the
            !! Dynamic Procedure. DON'T USE duidxjC again within this time step.
            !! Make the SGS call at the very end, just before the time
            !! advancement.
            
            if (this%isStratified) then
                call this%SGSmodel%getRHS_SGS_Scalar_WallM(this%dTdxC, this%dTdyC, this%dTdzE, &
                                                           this%T_rhs, this%wTh_surf           )
            end if 
        end if 
        
        ! Step 5: Time Step 
        if (this%step == 0) then
            this%uhat = this%uhat + this%dt*this%u_rhs 
            this%vhat = this%vhat + this%dt*this%v_rhs 
            this%what = this%what + this%dt*this%w_rhs 
            if (this%isStratified) then
                this%That = this%That + this%dt*this%T_rhs
            end if 
        else
            abf1 = (one + half*this%dtRat)*this%dt
            abf2 = -half*this%dtRat*this%dt
            this%uhat = this%uhat + abf1*this%u_rhs + abf2*this%u_Orhs
            this%vhat = this%vhat + abf1*this%v_rhs + abf2*this%v_Orhs
            this%what = this%what + abf1*this%w_rhs + abf2*this%w_Orhs
            if (this%isStratified) then
                this%That = this%That + abf1*this%T_rhs + abf2*this%T_Orhs
            end if 
        end if 
       
        
        ! Step 5: Dealias
        call this%spectC%dealias(this%uhat)
        call this%spectC%dealias(this%vhat)
        call this%spectE%dealias(this%what)
        if (this%isStratified) call this%spectC%dealias(this%That)
        if (this%UseDealiasFilterVert) then
            call this%ApplyCompactFilter()
        end if
       
        ! Step 6: Pressure projection
        if (useCompactFD) then
            call this%padepoiss%PressureProjection(this%uhat,this%vhat,this%what)
            if (mod(this%step,this%t_DivergenceCheck) == 0) then
                call this%padepoiss%DivergenceCheck(this%uhat, this%vhat, this%what, this%divergence)
            end if 
        else
            call this%poiss%PressureProjNP(this%uhat,this%vhat,this%what)
            if (mod(this%step,this%t_DivergenceCheck) == 0) then
                call this%poiss%DivergenceCheck(this%uhat, this%vhat, this%what, this%divergence)
            end if 
        end if 

        ! Step 7: Take it back to physical fields
        call this%spectC%ifft(this%uhat,this%u)
        call this%spectC%ifft(this%vhat,this%v)
        call this%spectE%ifft(this%what,this%w)
        if (this%isStratified) call this%spectC%ifft(this%That,this%T)
    
        ! STEP 8: Interpolate the cell center values of w
        call this%compute_and_bcast_surface_Mn()
        call this%interp_PrimitiveVars()

        ! STEP 9: Compute duidxjC 
        call this%compute_duidxj()
        if (this%isStratified) call this%compute_dTdxi() 
        
        ! STEP 10: Copy the RHS for using during next time step 
        this%u_Orhs = this%u_rhs
        this%v_Orhs = this%v_rhs
        this%w_Orhs = this%w_rhs
        if (this%isStratified) this%T_Orhs = this%T_rhs
        this%dtOld = this%dt


        ! STEP 11: Do logistical stuff
        if ((mod(this%step,this%tid_compStats)==0) .and. (this%tsim > this%tSimStartStats)) then
            call this%compute_stats()
        end if 

        if ((mod(this%step,this%tid_statsDump) == 0) .and. (this%tsim > this%tSimStartStats)) then
            call this%compute_stats()
            call this%dump_stats()
        end if 
        
        if (mod(this%step,this%t_restartDump) == 0) then
            call this%dumpRestartfile()
        end if
        
        if ((this%dumpPlanes) .and. (mod(this%step,this%t_planeDump) == 0) .and. &
                 (this%step .ge. this%t_start_planeDump) .and. (this%step .le. this%t_stop_planeDump)) then
            call this%dump_planes()
        end if 


        ! STEP 12: Update Time and BCs
        this%step = this%step + 1; this%tsim = this%tsim + this%dt
        if (this%isStratified) then
            this%Tsurf = this%Tsurf + this%dTsurf_dt*this%dt
        end if  

    end subroutine

    subroutine ApplyCompactFilter(this)
        class(igridWallM), intent(inout), target :: this
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
        call this%filzC%filter3(zbuff3,zbuff4,this%nxZ, this%nyZ)
        call transpose_z_to_y(zbuff4,this%what, this%sp_gpE)

        if (this%isStratified) then
            call transpose_y_to_z(this%That,zbuff1, this%sp_gpC)
            call this%filzC%filter3(zbuff1,zbuff2,this%nxZ, this%nyZ)
            call transpose_z_to_y(zbuff1,this%That, this%sp_gpC)
        end if

        nullify(zbuff1, zbuff2, zbuff3, zbuff4)
    end subroutine




    subroutine compute_duidxj(this)
        class(igridWallM), intent(inout), target :: this
        complex(rkind), dimension(:,:,:), pointer :: ctmpz1, ctmpz2
        complex(rkind), dimension(:,:,:), pointer :: ctmpz3!, ctmpz4
        complex(rkind), dimension(:,:,:), pointer :: ctmpy1, ctmpy2
        real(rkind),    dimension(:,:,:), pointer :: dudx, dudy, dudz
        real(rkind),    dimension(:,:,:), pointer :: dvdx, dvdy, dvdz
        real(rkind),    dimension(:,:,:), pointer :: dwdx, dwdy, dwdz
        real(rkind),    dimension(:,:,:), pointer :: dvdzC, dudzC
        real(rkind),    dimension(:,:,:), pointer :: dwdxC, dwdyC
        
        complex(rkind), dimension(:,:,:), pointer :: dudxH, dudyH, dudzH 
        complex(rkind), dimension(:,:,:), pointer :: dvdxH, dvdyH, dvdzH
        complex(rkind), dimension(:,:,:), pointer :: dwdxH, dwdyH, dwdzH

        complex(rkind), dimension(:,:)  , pointer :: dudz_dzby2, dvdz_dzby2

        dudx  => this%duidxjC(:,:,:,1); dudy  => this%duidxjC(:,:,:,2); dudzC => this%duidxjC(:,:,:,3); 
        dvdx  => this%duidxjC(:,:,:,4); dvdy  => this%duidxjC(:,:,:,5); dvdzC => this%duidxjC(:,:,:,6); 
        dwdxC => this%duidxjC(:,:,:,7); dwdyC => this%duidxjC(:,:,:,8); dwdz  => this%duidxjC(:,:,:,9); 

        dudxH => this%duidxjChat(:,:,:,1); dudyH => this%duidxjChat(:,:,:,2); dudzH => this%duidxjChat(:,:,:,3) 
        dvdxH => this%duidxjChat(:,:,:,4); dvdyH => this%duidxjChat(:,:,:,5); dvdzH => this%duidxjChat(:,:,:,6) 
        dwdxH => this%duidxjChat(:,:,:,7); dwdyH => this%duidxjChat(:,:,:,8); dwdzH => this%duidxjChat(:,:,:,9); 
       
        dwdx => this%duidxjE(:,:,:,1); dwdy => this%duidxjE(:,:,:,2);
        dudz => this%duidxjE(:,:,:,3); dvdz => this%duidxjE(:,:,:,4);

        ctmpz1 => this%cbuffzC(:,:,:,1); ctmpz2 => this%cbuffzE(:,:,:,1); 
        ctmpz3 => this%cbuffzC(:,:,:,2)!; ctmpz4 => this%cbuffzE(:,:,:,2)
        ctmpy1 => this%cbuffyC(:,:,:,1); ctmpy2 => this%cbuffyE(:,:,:,1)


        dudz_dzby2 => this%cbuffzE(:,:,1,2); dvdz_dzby2 => this%cbuffzE(:,:,1,2)


        call this%spectC%mTimes_ik1_oop(this%uhat,dudxH)
        call this%spectC%ifft(dudxH,dudx)

        call this%spectC%mTimes_ik2_oop(this%uhat,dudyH)
        call this%spectC%ifft(dudyH,dudy)

        call this%spectC%mTimes_ik1_oop(this%vhat,dvdxH)
        call this%spectC%ifft(dvdxH,dvdx)

        call this%spectC%mTimes_ik2_oop(this%vhat,dvdyH)
        call this%spectC%ifft(dvdyH,dvdy)
        
        call this%spectC%mTimes_ik1_oop(this%whatC,dwdxH)
        call this%spectC%ifft(dwdxH,dwdxC)

        call this%spectC%mTimes_ik2_oop(this%whatC,dwdyH)
        call this%spectC%ifft(dwdyH,dwdyC)

        call this%spectE%mTimes_ik1_oop(this%what,ctmpy2)
        call this%spectE%ifft(ctmpy2,dwdx)

        call this%spectE%mTimes_ik2_oop(this%what,ctmpy2)
        call this%spectE%ifft(ctmpy2,dwdy)
        
        call transpose_y_to_z(this%what,ctmpz2,this%sp_gpE)
        if (useCompactFD) then
            call this%derW%ddz_E2C(ctmpz2,ctmpz1,size(ctmpz2,1),size(ctmpz2,2))
        else
            call this%Ops%ddz_E2C(ctmpz2,ctmpz1)
        end if 
        call transpose_z_to_y(ctmpz1,dwdzH,this%sp_gpC)
        call this%spectC%ifft(dwdzH,dwdz)


        ! Compute dudz, dvdz 
        call transpose_y_to_z(this%uhat,ctmpz1,this%sp_gpC)
        if (useCompactFD) then
            call this%derEE%ddz_C2E(ctmpz1,ctmpz2,size(ctmpz1,1),size(ctmpz1,2))
        else    
            call this%Ops%ddz_C2E(ctmpz1,ctmpz2,topBC_u,.true.)
        end if     

        dudz_dzby2 = ctmpz1(:,:,1)/((this%dz/two)*log(this%dz/two/this%z0))
        call this%spectE%SurfaceFilter_ip(dudz_dzby2)
        dudz_dzby2 = (this%Umn/this%Uspmn) * dudz_dzby2

        ! Correct derivative at the z = dz (see Porte Agel, JFM (appendix))
        if (nrank == 0) then
            ctmpz2(1,1,2) = ctmpz2(1,1,2) + (0.08976d0/(kappa*this%dz))*real(this%nx,rkind)*real(this%ny,rkind)
        end if 


        ctmpz2(:,:,1) = (two*ctmpz2(:,:,2) - ctmpz2(:,:,3))
        call transpose_z_to_y(ctmpz2,ctmpy2,this%sp_gpE)
        call this%spectE%ifft(ctmpy2,dudz)
        if (useCompactFD) then
            call this%derOO%InterpZ_E2C(ctmpz2,ctmpz1,size(ctmpz2,1),size(ctmpz2,2))
        else
            call this%Ops%InterpZ_Edge2Cell(ctmpz2,ctmpz1)
        end if 
        ctmpz1(:,:,1) = dudz_dzby2 
        call transpose_z_to_y(ctmpz1,dudzH,this%sp_gpC)
        call this%spectC%ifft(dudzH,dudzC)

        call transpose_y_to_z(this%vhat,ctmpz1,this%sp_gpC)
        if (useCompactFD) then
            call this%derEE%ddz_C2E(ctmpz1,ctmpz2,size(ctmpz1,1),size(ctmpz1,2))
        else
            call this%Ops%ddz_C2E(ctmpz1,ctmpz2,topBC_v,.true.)
        end if 

        dvdz_dzby2 = dudz_dzby2 * this%Vmn/this%Umn
        ctmpz2(:,:,1) = two*dvdz_dzby2 - ctmpz2(:,:,2)
        call transpose_z_to_y(ctmpz2,ctmpy2,this%sp_gpE)
        call this%spectE%ifft(ctmpy2,dvdz)
        if (useCompactFD) then
            call this%derOO%InterpZ_E2C(ctmpz2,ctmpz1,size(ctmpz2,1),size(ctmpz2,2))
        else
            call this%Ops%InterpZ_Edge2Cell(ctmpz2,ctmpz1)
        end if 
        ctmpz1(:,:,1) = dvdz_dzby2
        call transpose_z_to_y(ctmpz1,dvdzH,this%sp_gpC)
        call this%spectC%ifft(dvdzH,dvdzC)

    end subroutine


    subroutine compute_dTdxi(this)
        class(igridWallM), intent(inout), target :: this
        complex(rkind), dimension(:,:,:), pointer :: ctmpz1, ctmpz2
        
        ctmpz1 => this%cbuffzC(:,:,:,1); ctmpz2 => this%cbuffzE(:,:,:,1); 
        
        call this%spectC%mtimes_ik1_oop(this%That,this%dTdxH)
        call this%spectC%ifft(this%dTdxH,this%dTdxC)

        call this%spectC%mtimes_ik2_oop(this%That,this%dTdyH)
        call this%spectC%ifft(this%dTdyH,this%dTdyC)
   
        call transpose_y_to_z(this%That, ctmpz1, this%sp_gpC)
        call this%OpsPP%ddz_C2E(ctmpz1,ctmpz2,.true.,.true.)
        ctmpz2(:,:,this%nz+1) = ctmpz2(:,:,this%nz)
        ctmpz2(:,:,1) = two*ctmpz2(:,:,2) - ctmpz2(:,:,3) 
        call transpose_z_to_y(ctmpz2, this%dTdzH, this%sp_gpE)
        call this%spectE%ifft(this%dTdzH,this%dTdzE)

    end subroutine
    
    subroutine debug(this)
        class(igridWallM), intent(inout), target :: this
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

        dwdx => this%duidxjE(:,:,:,1); dwdy => this%duidxjE(:,:,:,2);
        dudz => this%duidxjE(:,:,:,3); dvdz => this%duidxjE(:,:,:,4);

        rbuff => this%rbuffxC(:,:,:,1); cbuff => this%cbuffyC(:,:,:,1)
        dvdzH => this%duidxjChat(:,:,:,6) 

        !call this%spectC%ifft(this%v_rhs,rbuff)
        if (nrank == 0) then
            print*, "=------"
            print*, real(dvdzH(1,1,55:64))
        end if 
        !stop 
    end subroutine 
    
    subroutine compute_and_bcast_surface_Mn(this)
        use mpi
        use constants, only: four
        use kind_parameters, only: mpirkind
        class(igridWallM), intent(inout), target :: this
        real(rkind), dimension(:,:,:), pointer :: rbuff
        complex(rkind), dimension(:,:,:), pointer :: cbuff
        integer :: ierr

        rbuff => this%rbuffxC(:,:,:,1)
        cbuff => this%cbuffyC(:,:,:,1)
        rbuff = this%u*this%u
        rbuff = rbuff + this%v*this%v
        rbuff = sqrt(rbuff)
        call this%spectC%fft(rbuff,cbuff)

        if (nrank == 0) then
            this%Umn = real(this%uhat(1,1,1),rkind)*this%meanFact
            this%Vmn = real(this%vhat(1,1,1),rkind)*this%meanFact
            this%Uspmn = real(cbuff(1,1,1),rkind)*this%meanFact
            if (this%isStratified) this%Tmn = real(this%That(1,1,1),rkind)*this%meanFact
        end if
        call mpi_bcast(this%Umn,1,mpirkind,0,mpi_comm_world,ierr)
        call mpi_bcast(this%Vmn,1,mpirkind,0,mpi_comm_world,ierr)
        call mpi_bcast(this%Uspmn,1,mpirkind,0,mpi_comm_world,ierr)
        if (this%isStratified) call mpi_bcast(this%Tmn,1,mpirkind,0,mpi_comm_world,ierr)

        call this%getfilteredSpeedSqAtWall()
        if (this%isStratified) then
            call this%getSurfaceQuantities() 
        else
            this%ustar = this%Uspmn*kappa/(log(this%dz/two/this%z0))
            this%invObLength = zero
        end if 
    end subroutine
    
    subroutine getSurfaceQuantities(this)
        class(igridWallM), intent(inout) :: this
        integer :: idx
        integer, parameter :: itermax = 100 
        real(rkind) :: ustarNew, ustarDiff, dTheta, ustar
        real(rkind) :: a, b, c, PsiH, PsiM, wTh, z, u, Linv, g
        real(rkind), parameter :: beta_h = 7.8_rkind, beta_m = 4.8_rkind
       
        dTheta = this%Tsurf - this%Tmn
        z = this%dz/two ; ustarDiff = one; g = this%gravity_nd
        a=log(z/this%z0); b=beta_h*this%dz/two; c=beta_m*this%dz/two 
        PsiM = zero; PsiH = zero; idx = 0; ustar = one; u = this%Uspmn
        
        ! Inside the do loop all the used variables are on the stored on the stack
        ! After the while loop these variables are copied to their counterparts
        ! on the heap (variables part of the derived type)
        do while ( (ustarDiff > 1d-12) .and. (idx < itermax))
            ustarNew = u*kappa/(a - PsiM)
            wTh = dTheta*ustarNew*kappa/(a - PsiH) 
            Linv = -kappa*g*wTh/(ustarNew**3)
            PsiM = -c*Linv; PsiH = -b*Linv;
            ustarDiff = abs((ustarNew - ustar)/ustarNew)
            ustar = ustarNew; idx = idx + 1
        end do 
        this%ustar = ustar; this%invObLength = Linv; this%wTh_surf = wTh
    end subroutine

    subroutine readRestartFile(this, tid, rid)
        use decomp_2d_io
        use mpi
        use exits, only: message
        class(igridWallM), intent(inout) :: this
        integer, intent(in) :: tid, rid
        character(len=clen) :: tempname, fname
        integer :: ierr

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

        write(tempname,"(A7,A4,I2.2,A6,I6.6)") "RESTART", "_Run",rid, "_info.",tid
        fname = this%InputDir(:len_trim(this%InputDir))//"/"//trim(tempname)
        open(unit=10,file=fname,access='sequential',form='formatted')
        read (10, *)  this%tsim
        close(10)

        call mpi_barrier(mpi_comm_world, ierr)
        call message("================= RESTART FILE USED ======================")
        call message(0, "Simulation Time at restart:", this%tsim)
        call message("=================================== ======================")

    end subroutine

    subroutine dumpRestartFile(this)
        use decomp_2d_io
        use mpi
        use exits, only: message
        class(igridWallM), intent(in) :: this
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

        if (this%isStratified) then
            write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",this%runID, "_T.",this%step
            fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
            call decomp_2d_write_one(1,this%T,fname, this%gpE)
        end if 

        write(tempname,"(A7,A4,I2.2,A6,I6.6)") "RESTART", "_Run",this%runID, "_info.",this%step
        fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
        OPEN(UNIT=10, FILE=trim(fname))
        write(10,"(100g15.5)") this%tsim
        close(10)

        call mpi_barrier(mpi_comm_world, ierr)
        call message(1, "Just Dumped a RESTART file")

    end subroutine 


    subroutine dumpFullField(this,arr,label)
        use decomp_2d_io
        use mpi
        use exits, only: message
        class(igridWallM), intent(in) :: this
        character(len=clen) :: tempname, fname
        real(rkind), dimension(:,:,:), intent(in) :: arr
        character(len=4), intent(in) :: label

        write(tempname,"(A3,I2.2,A1,A4,A2,I6.6,A4)") "Run",this%runID, "_",label,"_t",this%step,".out"
        fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
        call decomp_2d_write_one(1,arr,fname)

    end subroutine


    !! STATISTICS !!

    subroutine init_stats( this)
        class(igridWallM), intent(inout), target :: this
        type(decomp_info), pointer  :: gpC

        gpC => this%gpC
        this%tidSUM = 0

        allocate(this%zStats2dump(this%nz,25))
        allocate(this%runningSum(this%nz,25))
        allocate(this%TemporalMnNOW(this%nz,25))

        ! mean velocities
        this%u_mean => this%zStats2dump(:,1);  this%v_mean  => this%zStats2dump(:,2);  this%w_mean => this%zStats2dump(:,3) 

        ! mean squared velocities
        this%uu_mean => this%zStats2dump(:,4); this%uv_mean => this%zStats2dump(:,5); this%uw_mean => this%zStats2dump(:,6)
                                               this%vv_mean => this%zStats2dump(:,7); this%vw_mean => this%zStats2dump(:,8) 
                                                                                      this%ww_mean => this%zStats2dump(:,9)

        ! SGS stresses
        this%tau11_mean => this%zStats2dump(:,10); this%tau12_mean => this%zStats2dump(:,11); this%tau13_mean => this%zStats2dump(:,12)
                                                   this%tau22_mean => this%zStats2dump(:,13); this%tau23_mean => this%zStats2dump(:,14) 
                                                                                              this%tau33_mean => this%zStats2dump(:,15)

        ! SGS dissipation
        this%sgsdissp => this%zStats2dump(:,16)

        ! velocity derivative products - for viscous dissipation
        this%viscdissp => this%zStats2dump(:,17)

        ! means of velocity derivatives
        this%S11_mean => this%zStats2dump(:,18); this%S12_mean => this%zStats2dump(:,19); this%S13_mean => this%zStats2dump(:,20)
                                                 this%S22_mean => this%zStats2dump(:,21); this%S23_mean => this%zStats2dump(:,22)
                                                                                          this%S33_mean => this%zStats2dump(:,23)

        ! SGS model coefficient
        this%sgscoeff_mean => this%zStats2dump(:,24)

        this%PhiM => this%zStats2dump(:,25)
        this%runningSum = zero
        nullify(gpC)
    end subroutine

    subroutine compute_stats(this)
        class(igridWallM), intent(inout), target :: this
        type(decomp_info), pointer :: gpC
        real(rkind), dimension(:,:,:), pointer :: rbuff1, rbuff2, rbuff3E, rbuff2E, rbuff3, rbuff4, rbuff5, rbuff5E, rbuff4E, rbuff6E, rbuff6

        rbuff1  => this%rbuffxC(:,:,:,1); rbuff2  => this%rbuffyC(:,:,:,1);
        rbuff2E => this%rbuffyE(:,:,:,1); rbuff3E => this%rbuffzE(:,:,:,1);
        rbuff3 => this%rbuffzC(:,:,:,1); rbuff4E => this%rbuffzE(:,:,:,2);
        rbuff4 => this%rbuffzC(:,:,:,2); rbuff5E => this%rbuffzE(:,:,:,3)
        rbuff5 => this%rbuffzC(:,:,:,3); rbuff6E => this%rbuffzE(:,:,:,4)
        rbuff6 => this%rbuffzC(:,:,:,4); 
        !rbuff7 => this%rbuffzC(:,:,:,5); 
        gpC => this%gpC

        this%tidSUM = this%tidSUM + 1

        ! Compute u - mean 
        call transpose_x_to_y(this%u,rbuff2,this%gpC)
        call transpose_y_to_z(rbuff2,rbuff3,this%gpC)
        call this%compute_z_mean(rbuff3, this%u_mean)
        if (this%normByustar) this%u_mean = this%u_mean/this%ustar

        ! Compute v - mean 
        call transpose_x_to_y(this%v,rbuff2,this%gpC)
        call transpose_y_to_z(rbuff2,rbuff4,this%gpC)
        call this%compute_z_mean(rbuff4, this%v_mean)
        if (this%normByustar)this%v_mean = this%v_mean/this%ustar

        ! Compute wC - mean 
        call transpose_x_to_y(this%wC,rbuff2,this%gpC)
        call transpose_y_to_z(rbuff2,rbuff5,this%gpC)
        call this%compute_z_mean(rbuff5, this%w_mean)
        if (this%normByustar)this%w_mean = this%w_mean/this%ustar

        ! take w from x -> z decomp
        call transpose_x_to_y(this%w,rbuff2E,this%gpE)
        call transpose_y_to_z(rbuff2E,rbuff5E,this%gpE)

        ! take uE from x -> z decomp
        call transpose_x_to_y(this%uE,rbuff2E,this%gpE)
        call transpose_y_to_z(rbuff2E,rbuff3E,this%gpE)

        ! take vE from x -> z decomp
        call transpose_x_to_y(this%vE,rbuff2E,this%gpE)
        call transpose_y_to_z(rbuff2E,rbuff4E,this%gpE)

        ! uu mean
        rbuff6 = rbuff3*rbuff3
        call this%compute_z_mean(rbuff6, this%uu_mean)
        if (this%normByustar)this%uu_mean = this%uu_mean/(this%ustar**2)

        ! uv mean
        rbuff6 = rbuff3*rbuff4
        call this%compute_z_mean(rbuff6, this%uv_mean)
        if (this%normByustar)this%uv_mean = this%uv_mean/(this%ustar**2)

        ! uw mean
        rbuff6E = rbuff3E*rbuff5E
        call this%OpsPP%InterpZ_Edge2Cell(rbuff6E,rbuff6)
        call this%compute_z_mean(rbuff6, this%uw_mean)
        if (this%normByustar)this%uw_mean = this%uw_mean/(this%ustar**2)

        ! vv mean 
        rbuff6 = rbuff4*rbuff4
        call this%compute_z_mean(rbuff6, this%vv_mean)
        if (this%normByustar)this%vv_mean = this%vv_mean/(this%ustar**2)

        ! vw mean 
        rbuff6E = rbuff4E*rbuff5E
        call this%OpsPP%InterpZ_Edge2Cell(rbuff6E,rbuff6)
        call this%compute_z_mean(rbuff6, this%vw_mean)
        if (this%normByustar)this%vw_mean = this%vw_mean/(this%ustar**2)

        ! ww mean 
        rbuff6 = rbuff5*rbuff5
        call this%compute_z_mean(rbuff6, this%ww_mean)
        if (this%normByustar)this%ww_mean = this%ww_mean/(this%ustar**2)

        if (this%useSGS) then
            ! tau_11
            call transpose_x_to_y(this%tauSGS_ij(:,:,:,1),rbuff2,this%gpC)
            call transpose_y_to_z(rbuff2,rbuff3,this%gpC)
            call this%compute_z_mean(rbuff3, this%tau11_mean)
            if (this%normByustar)this%tau11_mean = this%tau11_mean/(this%ustar**2)

            ! tau_12
            call transpose_x_to_y(this%tauSGS_ij(:,:,:,2),rbuff2,this%gpC)
            call transpose_y_to_z(rbuff2,rbuff3,this%gpC)
            call this%compute_z_mean(rbuff3, this%tau12_mean)
            if (this%normByustar)this%tau12_mean = this%tau12_mean/(this%ustar**2)

            ! tau_13
            call transpose_x_to_y(this%tau13,rbuff2E,this%gpE)
            call transpose_y_to_z(rbuff2E,rbuff3E,this%gpE)
            rbuff3E(:,:,1) = -(this%ustar**2)
            call this%OpsPP%InterpZ_Edge2Cell(rbuff3E,rbuff3)
            call this%compute_z_mean(rbuff3, this%tau13_mean)
            if (this%normByustar)this%tau13_mean = this%tau13_mean/(this%ustar**2)

            ! tau_22
            call transpose_x_to_y(this%tauSGS_ij(:,:,:,4),rbuff2,this%gpC)
            call transpose_y_to_z(rbuff2,rbuff3,this%gpC)
            call this%compute_z_mean(rbuff3, this%tau22_mean)
            if (this%normByustar)this%tau22_mean = this%tau22_mean/(this%ustar**2)

            ! tau_23
            call transpose_x_to_y(this%tau23,rbuff2E,this%gpE)
            call transpose_y_to_z(rbuff2E,rbuff3E,this%gpE)
            rbuff3E(:,:,1) = -(this%ustar**2)*this%Vmn/this%Umn
            call this%OpsPP%InterpZ_Edge2Cell(rbuff3E,rbuff3)
            call this%compute_z_mean(rbuff3, this%tau23_mean)
            if (this%normByustar)this%tau23_mean = this%tau23_mean/(this%ustar**2)

            ! tau_33
            call transpose_x_to_y(this%tauSGS_ij(:,:,:,6),rbuff2,this%gpC)
            call transpose_y_to_z(rbuff2,rbuff3,this%gpC)
            call this%compute_z_mean(rbuff3, this%tau33_mean)
            if (this%normByustar)this%tau33_mean = this%tau33_mean/(this%ustar**2)


            ! sgs dissipation
            rbuff1 = this%tauSGS_ij(:,:,:,1)*this%tauSGS_ij(:,:,:,1) + &
                     this%tauSGS_ij(:,:,:,2)*this%tauSGS_ij(:,:,:,2) + &
                     this%tauSGS_ij(:,:,:,3)*this%tauSGS_ij(:,:,:,3)
            rbuff1 = rbuff1 + two*(this%tauSGS_ij(:,:,:,4)*this%tauSGS_ij(:,:,:,4) + &
                                   this%tauSGS_ij(:,:,:,5)*this%tauSGS_ij(:,:,:,5) + &
                                   this%tauSGS_ij(:,:,:,6)*this%tauSGS_ij(:,:,:,6) )
            rbuff1 = rbuff1/(this%nu_SGS + 1.0d-14)         ! note: factor of half is in dump_stats

            call transpose_x_to_y(rbuff1,rbuff2,this%gpC)
            call transpose_y_to_z(rbuff2,rbuff3,this%gpC)
            call this%compute_z_mean(rbuff3, this%sgsdissp)

            ! viscous dissipation- *****????? Is rbuff1 contaminated after transpose_x_to_y? *****?????
            rbuff1 = rbuff1/(this%nu_SGS + 1.0d-14)        ! note: factor of fourth is in dump_stats
            call transpose_x_to_y(rbuff1,rbuff2,this%gpC)
            call transpose_y_to_z(rbuff2,rbuff3,this%gpC)
            call this%compute_z_mean(rbuff3, this%viscdissp)

            ! note: factor of half in all S_** is in dump_stats
            ! S_11
            rbuff1 = this%tauSGS_ij(:,:,:,1)/(this%nu_SGS + 1.0d-14)
            call transpose_x_to_y(rbuff1,rbuff2,this%gpC)
            call transpose_y_to_z(rbuff2,rbuff3,this%gpC)
            call this%compute_z_mean(rbuff3, this%S11_mean)

            ! S_12
            rbuff1 = this%tauSGS_ij(:,:,:,2)/(this%nu_SGS + 1.0d-14)
            call transpose_x_to_y(rbuff1,rbuff2,this%gpC)
            call transpose_y_to_z(rbuff2,rbuff3,this%gpC)
            call this%compute_z_mean(rbuff3, this%S12_mean)

            ! S_13
            rbuff1 = this%tauSGS_ij(:,:,:,3)/(this%nu_SGS + 1.0d-14)
            call transpose_x_to_y(rbuff1,rbuff2,this%gpC)
            call transpose_y_to_z(rbuff2,rbuff3,this%gpC)
            call this%compute_z_mean(rbuff3, this%S13_mean)

            ! S_22
            rbuff1 = this%tauSGS_ij(:,:,:,4)/(this%nu_SGS + 1.0d-14)
            call transpose_x_to_y(rbuff1,rbuff2,this%gpC)
            call transpose_y_to_z(rbuff2,rbuff3,this%gpC)
            call this%compute_z_mean(rbuff3, this%S22_mean)

            ! S_23
            rbuff1 = this%tauSGS_ij(:,:,:,5)/(this%nu_SGS + 1.0d-14)
            call transpose_x_to_y(rbuff1,rbuff2,this%gpC)
            call transpose_y_to_z(rbuff2,rbuff3,this%gpC)
            call this%compute_z_mean(rbuff3, this%S23_mean)

            ! S_33
            rbuff1 = this%tauSGS_ij(:,:,:,6)/(this%nu_SGS + 1.0d-14)
            call transpose_x_to_y(rbuff1,rbuff2,this%gpC)
            call transpose_y_to_z(rbuff2,rbuff3,this%gpC)
            call this%compute_z_mean(rbuff3, this%S33_mean)

            ! sgs coefficient
            call transpose_x_to_y(this%c_SGS,rbuff2,this%gpC)
            call transpose_y_to_z(rbuff2,rbuff3,this%gpC)
            !call this%compute_z_mean(rbuff3, this%sgscoeff_mean)    ! -- averaging not needed
            this%sgscoeff_mean(:) = rbuff3(1,1,:)
        
        end if

        rbuff1 = this%duidxjC(:,:,:,3)*this%mesh(:,:,:,3)
        call transpose_x_to_y(rbuff1,rbuff2,this%gpC)
        call transpose_y_to_z(rbuff2,rbuff3,this%gpC)
        call this%compute_z_mean(rbuff3, this%PhiM)
        this%PhiM = this%PhiM*kappa/this%ustar

        this%runningSum = this%runningSum + this%zStats2dump

        !write(*,*) 'In stats'
        !write(*,*) 'umean', maxval(this%u_mean), minval(this%u_mean)
        !write(*,*) 'vmean', maxval(this%v_mean), minval(this%v_mean)
        !write(*,*) 'wmean', maxval(this%w_mean), minval(this%w_mean)

    end subroutine 

    subroutine dump_stats(this)
        use basic_io, only: write_2d_ascii, write_2D_binary
        use exits, only: message
        use kind_parameters, only: clen
        use mpi
        class(igridWallM), intent(inout), target :: this
        character(len=clen) :: fname
        character(len=clen) :: tempname
        integer :: tid

        this%TemporalMnNOW = this%runningSum/real(this%tidSUM,rkind)
        tid = this%step

        ! compute (u_i'u_j')
        this%TemporalMnNOW(:,4) = this%TemporalMnNOW(:,4) - this%TemporalMnNOW(:,1)*this%TemporalMnNOW(:,1)
        this%TemporalMnNOW(:,5) = this%TemporalMnNOW(:,5) - this%TemporalMnNOW(:,1)*this%TemporalMnNOW(:,2)
        this%TemporalMnNOW(:,6) = this%TemporalMnNOW(:,6) - this%TemporalMnNOW(:,1)*this%TemporalMnNOW(:,3)
        this%TemporalMnNOW(:,7) = this%TemporalMnNOW(:,7) - this%TemporalMnNOW(:,2)*this%TemporalMnNOW(:,2)
        this%TemporalMnNOW(:,8) = this%TemporalMnNOW(:,8) - this%TemporalMnNOW(:,2)*this%TemporalMnNOW(:,3)
        this%TemporalMnNOW(:,9) = this%TemporalMnNOW(:,9) - this%TemporalMnNOW(:,3)*this%TemporalMnNOW(:,3)

        ! compute sgs dissipation
        this%TemporalMnNOW(:,16) = half*this%TemporalMnNOW(:,16)

        ! compute viscous dissipation
        this%TemporalMnNOW(:,17) = this%TemporalMnNOW(:,17) - (                        &
                                   this%TemporalMnNOW(:,18)*this%TemporalMnNOW(:,18) + &
                                   this%TemporalMnNOW(:,21)*this%TemporalMnNOW(:,21) + &
                                   this%TemporalMnNOW(:,23)*this%TemporalMnNOW(:,23) + &
                              two*(this%TemporalMnNOW(:,19)*this%TemporalMnNOW(:,19) + &
                                   this%TemporalMnNOW(:,20)*this%TemporalMnNOW(:,20) + & 
                                   this%TemporalMnNOW(:,22)*this%TemporalMnNOW(:,22)))
        this%TemporalMnNOW(:,17) = half*this%TemporalMnNOW(:,17)/this%Re     ! note: this is actually 2/Re*(..)/4

        if (nrank == 0) then
            write(tempname,"(A3,I2.2,A2,I6.6,A4)") "Run", this%RunID,"_t",tid,".stt"
            fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
            call write_2d_ascii(this%TemporalMnNOW,fname)
            !call write_2D_binary(TemporalMnNOW,fname)
        end if
        call message(1, "Just dumped a .stt file")
        call message(2, "Number ot tsteps averaged:",this%tidSUM)

    end subroutine

    subroutine compute_z_fluct(this,fin)
        use reductions, only: P_SUM
        class(igridWallM), intent(in), target :: this
        real(rkind), dimension(:,:,:), intent(inout) :: fin
        integer :: k
        real(rkind) :: fmean

        do k = 1,size(fin,3)
            fmean = P_SUM(sum(fin(:,:,k)))/(real(this%nx,rkind)*real(this%ny,rkind))
            fin(:,:,k) = fin(:,:,k) - fmean
        end do 

    end subroutine

    subroutine compute_z_mean(this, arr_in, vec_out)
        use reductions, only: P_SUM
        class(igridWallM), intent(in), target :: this
        real(rkind), dimension(:,:,:), intent(in) :: arr_in
        real(rkind), dimension(:), intent(out) :: vec_out
        integer :: k

        do k = 1,size(arr_in,3)
            vec_out(k) = P_SUM(sum(arr_in(:,:,k)))/(real(this%nx,rkind)*real(this%ny,rkind))
        end do 

    end subroutine

    subroutine finalize_stats(this)
        class(igridWallM), intent(inout) :: this
        nullify(this%u_mean, this%v_mean, this%w_mean, this%uu_mean, this%uv_mean, this%uw_mean, this%vv_mean, this%vw_mean, this%ww_mean)
        nullify(this%tau11_mean, this%tau12_mean, this%tau13_mean, this%tau22_mean, this%tau23_mean, this%tau33_mean)
        nullify(this%S11_mean, this%S12_mean, this%S13_mean, this%S22_mean, this%S23_mean, this%S33_mean)
        nullify(this%sgsdissp, this%viscdissp, this%sgscoeff_mean)
        deallocate(this%zStats2dump, this%runningSum, this%TemporalMnNOW)
    end subroutine 

    subroutine dump_planes(this)
        use decomp_2d_io
        class(igridWallM), intent(in) :: this
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
                call decomp_2d_write_plane(1,this%u,dirid, pid, fname)

                write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_x",pid,".plv"
                fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                call decomp_2d_write_plane(1,this%v,dirid, pid, fname)

                write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_x",pid,".plw"
                fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                call decomp_2d_write_plane(1,this%wC,dirid, pid, fname)
            end do 
        end if 
            
            
        if (allocated(this%yplanes)) then
            nyplanes = size(this%yplanes)
            dirid = 2
            do idx = 1,nyplanes
                pid = this%yplanes(idx)
                write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_y",pid,".plu"
                fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                call decomp_2d_write_plane(1,this%u,dirid, pid, fname)

                write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_y",pid,".plv"
                fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                call decomp_2d_write_plane(1,this%v,dirid, pid, fname)

                write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_y",pid,".plw"
                fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                call decomp_2d_write_plane(1,this%wC,dirid, pid, fname)
            end do 
        end if 
        
        
        if (allocated(this%zplanes)) then
            nzplanes = size(this%zplanes)
            dirid = 3
            do idx = 1,nzplanes
                pid = this%zplanes(idx)
                write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_z",pid,".plu"
                fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                call decomp_2d_write_plane(1,this%u,dirid, pid, fname)

                write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_z",pid,".plv"
                fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                call decomp_2d_write_plane(1,this%v,dirid, pid, fname)

                write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_z",pid,".plw"
                fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                call decomp_2d_write_plane(1,this%wC,dirid, pid, fname)
            end do 
        end if 
        call message(1, "Dumped Planes.")        
    end subroutine 



    subroutine getfilteredSpeedSqAtWall(this)
        class(igridWallM), intent(inout), target :: this

        real(rkind), dimension(:,:,:), pointer :: rbuffx1, rbuffx2
        complex(rkind), dimension(:,:,:), pointer :: cbuffy, tauWallH

        cbuffy => this%cbuffyC(:,:,:,1); tauWallH => this%cbuffzC(:,:,:,1)     
        rbuffx1 => this%filteredSpeedSq; rbuffx2 => this%rbuffxC(:,:,:,1)

        call transpose_y_to_z(this%uhat,tauWallH,this%sp_gpC)
        call this%spectC%SurfaceFilter_ip(tauWallH(:,:,1))
        call transpose_z_to_y(tauWallH,cbuffy, this%sp_gpC)
        call this%spectC%ifft(cbuffy,rbuffx1)

        call transpose_y_to_z(this%vhat,tauWallH,this%sp_gpC)
        call this%spectC%SurfaceFilter_ip(tauWallH(:,:,1))
        call transpose_z_to_y(tauWallH,cbuffy, this%sp_gpC)
        call this%spectC%ifft(cbuffy,rbuffx2)

        rbuffx1 = rbuffx1*rbuffx1
        rbuffx2 = rbuffx2*rbuffx2
        rbuffx1 = rbuffx1 + rbuffx2

    end subroutine  



end module 
