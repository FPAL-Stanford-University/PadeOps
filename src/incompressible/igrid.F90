module IncompressibleGridNP
    use kind_parameters, only: rkind, clen
    use constants, only: imi, zero,one,two,three,half,fourth 
    use GridMod, only: grid
    use gridtools, only: alloc_buffs, destroy_buffs
    use igrid_hooks, only: meshgen, initfields_stagg, getforcing, set_planes_io
    use decomp_2d
    use StaggOpsMod, only: staggOps  
    use exits, only: GracefulExit, message
    use spectralMod, only: spectral  
    use PoissonMod, only: poisson
    use mpi 
    use reductions, only: p_maxval
    use timer, only: tic, toc
    use PadePoissonMod, only: Padepoisson 
    use cd06staggstuff, only: cd06stagg
    use cf90stuff, only: cf90  
    use sgsmod, only: sgs
    use numerics, only: useCompactFD, AdvectionForm
    use wallmodelMod, only: wallmodel

    implicit none

    private
    public :: igrid 


    complex(rkind), parameter :: zeroC = zero + imi*zero 
   integer, parameter :: no_slip = 1, slip = 2


    ! Allow non-zero value (isEven) 
    logical :: topBC_u = .true.  , topBC_v = .true. 
    logical :: botBC_u = .false. , botBC_v = .false. 
    logical, parameter :: topBC_w = .false. , botBC_w = .false. 
    integer :: ierr 

    type, extends(grid) :: igrid
        
        character(clen) :: inputDir

        type(decomp_info), allocatable :: gpC, gpE
        type(decomp_info), pointer :: Sp_gpC, Sp_gpE
        type(spectral), allocatable :: spectE, spectC
        type(staggOps), allocatable :: Ops
        type(sgs), allocatable :: SGSmodel
        type(wallmodel), allocatable :: moengWall

        real(rkind), dimension(:,:,:,:), allocatable :: PfieldsC
        real(rkind), dimension(:,:,:,:), allocatable :: PfieldsE

        complex(rkind), dimension(:,:,:,:), allocatable :: SfieldsC
        complex(rkind), dimension(:,:,:,:), allocatable :: SfieldsE


        type(poisson), allocatable :: poiss
        type(padepoisson), allocatable :: padepoiss
        real(rkind), dimension(:,:,:), allocatable :: divergence

        real(rkind), dimension(:,:,:), pointer :: u, v, wC, w
        real(rkind), dimension(:,:,:), pointer :: ox,oy,oz
        complex(rkind), dimension(:,:,:), pointer :: uhat, vhat, whatC, what
        complex(rkind), dimension(:,:,:), pointer :: oxhat, oyhat, ozhat

        type(cd06stagg), allocatable :: derU, derV, derW, derWW
        type(cf90),      allocatable :: filzE, filzC

        real(rkind), dimension(:,:,:,:), allocatable, public :: rbuffxC, rbuffyC, rbuffzC
        real(rkind), dimension(:,:,:,:), allocatable :: rbuffxE, rbuffyE, rbuffzE
        
        complex(rkind), dimension(:,:,:,:), allocatable :: cbuffyC, cbuffzC
        complex(rkind), dimension(:,:,:,:), allocatable :: cbuffyE, cbuffzE

        complex(rkind), dimension(:,:,:,:), allocatable :: rhsC, rhsE, OrhsC, OrhsE 
        real(rkind), dimension(:,:,:,:), allocatable :: duidxj 
        complex(rkind), dimension(:,:,:,:), allocatable :: duidxjhat
        complex(rkind), dimension(:,:,:), pointer:: u_rhs, v_rhs, wC_rhs, w_rhs 
        complex(rkind), dimension(:,:,:), pointer:: u_Orhs, v_Orhs, w_Orhs
            
        real(rkind) :: Re, Gx, Gy, Gz, dtby2
        complex(rkind), dimension(:,:,:), allocatable :: GxHat 
        real(rkind) :: Ro = 1.d5

        integer :: nxZ, nyZ
        integer :: tid_statsDump
        real(rkind) :: time_startDumping 
        
        integer :: runID
        logical :: useCoriolis = .true. 
        logical :: useExtraForcing = .false.
        logical :: useWallModelTop = .false., useWallModelBot = .false. 
        logical :: isInviscid = .false. 
        logical :: useSGS = .false. 
        logical :: useVerticalFilter = .true. 
        logical :: useDynamicProcedure 

        complex(rkind), dimension(:,:,:), allocatable :: dPf_dxhat

        real(rkind) :: max_nuSGS

        ! Statistics to compute 
        real(rkind), dimension(:,:), allocatable :: zStats2dump, runningSum, TemporalMnNOW
        real(rkind), dimension(:), pointer :: u_mean, v_mean, w_mean, uu_mean, uv_mean, uw_mean, vv_mean, vw_mean, ww_mean
        real(rkind), dimension(:), pointer :: tau11_mean, tau12_mean, tau13_mean, tau22_mean, tau23_mean, tau33_mean
        real(rkind), dimension(:), pointer :: S11_mean, S12_mean, S13_mean, S22_mean, S23_mean, S33_mean
        real(rkind), dimension(:), pointer :: viscdissp, sgsdissp, sgscoeff_mean
        integer :: tidSUM


        ! Pointers linked to SGS stuff
        real(rkind), dimension(:,:,:,:), pointer :: tauSGS_ij
        real(rkind), dimension(:,:,:)  , pointer :: nu_SGS
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
            !procedure :: laplacian
            !procedure :: gradient
            procedure :: AdamsBashforth
            procedure :: getMaxKE
            procedure, private :: interp_wHat_to_wHatC
            procedure :: compute_vorticity
            procedure, private :: compute_duidxj
            procedure, private :: addNonLinearTerm_Rot
            procedure, private :: addNonLinearTerm_SkewSymm 
            procedure, private :: addCoriolisTerm
            procedure, private :: addViscousTerm 
            procedure, private :: addExtraForcingTerm 
            procedure, private :: ApplyCompactFilter
            procedure          :: dumpRestartFile
            procedure, private :: readRestartFile
            procedure, private :: compute_z_mean 
            procedure, private :: compute_z_fluct
            procedure          :: dump_stats
            procedure          :: compute_stats 
            procedure          :: finalize_stats
            procedure          :: dump_planes 
    end type

contains 

    subroutine init(this,inputfile)
        class(igrid), intent(inout), target :: this        
        character(len=clen), intent(in) :: inputfile 
        integer :: nx, ny, nz
        character(len=clen) :: outputdir
        character(len=clen) :: inputdir
        logical :: SkewSymm = .false. 
        logical :: periodicx = .true. 
        logical :: periodicy = .true. 
        logical :: periodicz = .true.
        character(len=clen) :: derivative_x = "four"  
        character(len=clen) :: derivative_y = "four" 
        character(len=clen) :: derivative_z = "four"
        character(len=clen) :: filter_x = "2/3rd"  
        character(len=clen) :: filter_y = "2/3rd" 
        character(len=clen) :: filter_z = "2/3rd"
        integer :: prow = 0, pcol = 0 
        logical :: useSGS = .false., useDynamicProcedure = .false.
        logical :: useSGSclipping = .true.  
        integer :: ioUnit
        integer :: nsteps = -1
        real(rkind) :: dt = -one
        real(rkind) :: tstop = one
        real(rkind) :: CFL = -one
        real(rkind) :: Re = 800.00_rkind
        real(rkind) :: u_g = 1._rkind
        integer :: runID = 0
        integer :: t_dataDump = 99999
        integer :: tid_statsDump = 1000
        real(rkind) :: time_startDumping = 100000._rkind
        integer :: t_restartDump = 99999
        logical :: ViscConsrv = .TRUE. 
        real(rkind) :: Pr = 0.7_rkind 
        integer :: topWall = slip
        integer :: botWall = no_slip
        logical :: useCoriolis = .true. 
        logical :: useExtraForcing = .false.
        real(rkind) :: dpFdx = zero
        integer :: restartFile_TID = 1
        integer :: restartFile_RID = 1
        logical :: useRestartFile = .false. 
        logical :: useWallModelTop = .false., useWallModelBot = .false.
        logical :: isInviscid = .false., useVerticalFilter = .true.  
        integer :: SGSModelID = 1
        integer :: vscheme = 0, advForm = 1

        namelist /INPUT/       nx, ny, nz, tstop, dt, CFL, nsteps, &
                                              inputdir, outputdir, &
                                  periodicx, periodicy, periodicz, &
                                                       prow, pcol, &
                                        t_restartDump, t_dataDump
        namelist /IINPUT/  Re, useSGS, useDynamicProcedure,runID, Pr, useCoriolis, & 
                                tid_statsDump, useExtraForcing, useSGSclipping, &
                                time_startDumping, topWall, botWall, &
                                useRestartFile, restartFile_TID, restartFile_RID, &
                                isInviscid, &
                                useVerticalFilter, SGSModelID, vscheme, advForm 

        ! STEP 1: READ INPUT 
        ioUnit = 11
        open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
        read(unit=ioUnit, NML=INPUT)
        read(unit=ioUnit, NML=IINPUT)
        close(ioUnit)
      
        this%nx = nx
        this%ny = ny
        this%nz = nz

        this%SkewSymm = SkewSymm 
        this%ViscConsrv = ViscConsrv
        this%dt = dt
        this%dtby2 = dt/two 

        this%Re = Re
        this%outputdir = outputdir 
        this%inputdir = inputdir 
        
        this%periodicx = periodicx
        this%periodicy = periodicy
        this%periodicz = periodicz

        this%derivative_x = derivative_x    
        this%derivative_y = derivative_y    
        this%derivative_z = derivative_z  

        this%filter_x = filter_x    
        this%filter_y = filter_y    
        this%filter_z = filter_z  

        this%runID = runID
        this%tstop = tstop 
        this%t_dataDump = t_dataDump

        this%t_restartDump = t_restartDump
        this%tid_statsDump = tid_statsDump
        this%time_startDumping = time_startDumping 
        this%useCoriolis = useCoriolis 
        this%useExtraForcing = useExtraForcing

        this%useWallModelTop = useWallModelTop
        this%useWallModelBot = useWallModelBot
        this%isInviscid = isInviscid

        this%useSGS = useSGS
        this%UseDynamicProcedure = useDynamicProcedure

        this%useVerticalFilter = useVerticalFilter

        select case (vscheme)
        case(0)
            useCompactFD = .true. 
        case(1)
            useCompactFD = .false.
        case default 
            call GracefulExit("Invalid choice for VSCHEME. Only options &
                & available are 0 (6th Order) and 1 (2nd order)",3214)
        end select

        AdvectionForm = advForm 

        if (.not. periodicx) then
            call GracefulExit("Currently only Periodic BC is supported in x direction",102)
        end if 

        if (.not. periodicy) then
            call GracefulExit("Currently only Periodic BC is supported in y direction",102)
        end if 
        
        if (this%periodicz) then
            call GracefulExit("For all periodic direction problems, use HIT_GRID instead of IGRID",102)
        end if 

        
        ! STEP 2: ALLOCATE DECOMPOSITIONS
        allocate(this%gpC)
        allocate(this%gpE)

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
            call message("WARNING: No Top BCs provided. Using defaults found in igrid.F90")
        end if 
        
        if (botWall == slip) then
            botBC_u = .true.; botBC_v = .true.
            call message(1, "BotWall BC set to: SLIP")
        elseif (botWall == no_slip) then
            botBC_u = .false.; botBC_v = .false.
            call message(1, "BotWall BC set to: NO_SLIP")
        elseif (botWall == 3) then
            botBC_u = .true.; botBC_v = .true.
            useWallModelBot = .true. 
            call message(1, "BotWall BC set to: WALL MODEL (Moeng)")
        else
            call message("WARNING: No Bottom BCs provided. Using defaults found in igrid.F90")
        end if 
        
        
        ! STEP 3: GENERATE MESH (CELL CENTERED) 
        if ( allocated(this%mesh) ) deallocate(this%mesh) 
        allocate(this%mesh(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3),3))
        call meshgen(this%gpC, this%dx, this%dy, &
            this%dz, this%mesh) ! <-- this procedure is part of user defined HOOKS
        
        ! STEP 4: ALLOCATE/INITIALIZE THE SPECTRAL DERIVED TYPES
        allocate(this%spectC)
        call this%spectC%init("x", nx, ny, nz, this%dx, this%dy, this%dz, &
                this%derivative_x, this%filter_x, 2 , .false.)
        allocate(this%spectE)
        call this%spectE%init("x", nx, ny, nz+1, this%dx, this%dy, this%dz, &
                this%derivative_x, this%filter_x, 2 , .false.)
        this%sp_gpC => this%spectC%spectdecomp
        this%sp_gpE => this%spectE%spectdecomp


        ! STEP 5: ALLOCATE/INITIALIZE THE OPERATORS DERIVED TYPE
        if (useCompactFD) then
            allocate(this%derU, this%derV, this%derW, this%derWW)
            call this%derU%init (nz  , this%dz, topBC_u, botBC_u, .false., useWallModelBot) 
            call this%derV%init (nz  , this%dz, topBC_v, botBC_v, .false., useWallModelBot)  
            call this%derW%init (nz  , this%dz, topBC_w, botBC_w)
            call this%derWW%init(nz , this%dz, .true., .true.)
        else
            allocate(this%Ops)
            call this%Ops%init(this%gpC,this%gpE,0,this%dx,this%dy,this%dz,this%spectC%spectdecomp, &
                        this%spectE%spectdecomp, .false., useWallModelBot)
        end if 
        
        if (this%useVerticalFilter) then
            allocate(this%filzC, this%filzE)
            ierr = this%filzC%init(nz  , .false.)
            ierr = this%filzE%init(nz+1, .false.)
        end if 

        ! STEP 6: ALLOCATE MEMORY FOR FIELD ARRAYS
        allocate(this%PfieldsC(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3),6))
        allocate(this%divergence(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)))
        allocate(this%duidxj(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3),9))
        allocate(this%PfieldsE(this%gpE%xsz(1),this%gpE%xsz(2),this%gpE%xsz(3),1))
        call this%spectC%alloc_r2c_out(this%SfieldsC,6)
        call this%spectC%alloc_r2c_out(this%duidxjhat,9)
        call this%spectC%alloc_r2c_out(this%rhsC,3)
        call this%spectC%alloc_r2c_out(this%OrhsC,2)
        call this%spectE%alloc_r2c_out(this%rhsE,1)
        call this%spectE%alloc_r2c_out(this%OrhsE,1)

        
        call this%spectE%alloc_r2c_out(this%SfieldsE,1)
        
        this%u => this%PfieldsC(:,:,:,1) 
        this%v => this%PfieldsC(:,:,:,2) 
        this%wC => this%PfieldsC(:,:,:,3) 
        this%w => this%PfieldsE(:,:,:,1) 
        
        this%uhat => this%SfieldsC(:,:,:,1)
        this%vhat => this%SfieldsC(:,:,:,2)
        this%whatC => this%SfieldsC(:,:,:,3)
        this%what => this%SfieldsE(:,:,:,1)

        this%ox => this%PfieldsC(:,:,:,4)
        this%oy => this%PfieldsC(:,:,:,5)
        this%oz => this%PfieldsC(:,:,:,6)

        this%oxhat => this%SfieldsC(:,:,:,4)
        this%oyhat => this%SfieldsC(:,:,:,5)
        this%ozhat => this%SfieldsC(:,:,:,6)

        this%u_rhs => this%rhsC(:,:,:,1)
        this%v_rhs => this%rhsC(:,:,:,2)
        this%wC_rhs => this%rhsC(:,:,:,3)
        this%w_rhs => this%rhsE(:,:,:,1)

        this%u_Orhs => this%OrhsC(:,:,:,1)
        this%v_Orhs => this%OrhsC(:,:,:,2)
        this%w_Orhs => this%OrhsE(:,:,:,1)

        allocate(this%cbuffyC(this%sp_gpC%ysz(1),this%sp_gpC%ysz(2),this%sp_gpC%ysz(3),2))
        allocate(this%cbuffyE(this%sp_gpE%ysz(1),this%sp_gpE%ysz(2),this%sp_gpE%ysz(3),2))
        
        allocate(this%cbuffzC(this%sp_gpC%zsz(1),this%sp_gpC%zsz(2),this%sp_gpC%zsz(3),2))
        allocate(this%cbuffzE(this%sp_gpE%zsz(1),this%sp_gpE%zsz(2),this%sp_gpE%zsz(3),2))

        allocate(this%rbuffxC(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3),2))
        allocate(this%rbuffyC(this%gpC%ysz(1),this%gpC%ysz(2),this%gpC%ysz(3),2))
        allocate(this%rbuffzC(this%gpC%zsz(1),this%gpC%zsz(2),this%gpC%zsz(3),4))

        allocate(this%rbuffxE(this%gpE%xsz(1),this%gpE%xsz(2),this%gpE%xsz(3),2))
        allocate(this%rbuffyE(this%gpE%ysz(1),this%gpE%ysz(2),this%gpE%ysz(3),2))
        allocate(this%rbuffzE(this%gpE%zsz(1),this%gpE%zsz(2),this%gpE%zsz(3),2))
        this%nxZ = size(this%cbuffzE,1)
        this%nyZ = size(this%cbuffzE,2)


        ! STEP 6: ALLOCATE/INITIALIZE THE POISSON DERIVED TYPE
        if (useCompactFD) then
            allocate(this%padepoiss)
            call this%padepoiss%init(this%dx, this%dy, this%dz, this%spectC, this%spectE, this%derW) 
        else    
            allocate(this%poiss)
            call this%poiss%init(this%spectC,.false.,this%dx,this%dy,this%dz,this%Ops,this%spectE)  
        end if 
               
 
        ! STEP 7: INITIALIZE THE FIELDS
        if (useRestartFile) then
            call this%readRestartFile(restartfile_TID, restartfile_RID)
            this%step = restartfile_TID
        else 
            call initfields_stagg(this%gpC, this%gpE, this%dx, this%dy, this%dz, &
                inputfile, this%mesh, this%PfieldsC, this%PfieldsE, u_g, this%Ro)! <-- this procedure is part of user defined HOOKS
            this%step = 0
            this%tsim = zero
            call this%dumpRestartfile()
        end if 
        
        this%Gx = u_g
        this%Gy = zero
        this%Gz = zero
        

        call this%spectC%fft(this%u,this%uhat)   
        call this%spectC%fft(this%v,this%vhat)   
        call this%spectE%fft(this%w,this%what)   
      
        ! Dealias and filter before projection
        call this%spectC%dealias(this%uhat)
        call this%spectC%dealias(this%vhat)
        call this%spectE%dealias(this%what)
        if (this%useVerticalFilter) then
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
        
        ! STEP 8: Interpolate the cell center values of w
        call this%interp_wHat_to_wHatC()
        call message(1,"Max KE:",P_MAXVAL(this%getMaxKE()))
     
        ! STEP 9: Compute Vorticity
        if (AdvectionForm == 1) then 
            call this%compute_Vorticity()
        end if 
        if ((AdvectionForm == 2) .or. (this%useSGS)) then
            call this%compute_duidxj()
        end if 
        

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
            call getForcing(inputfile, dpFdx)
            call message(0," Turning on aditional forcing")
            call message(1," dP_dx = ", dpFdx)
            call this%spectC%alloc_r2c_out(this%dpF_dxhat)
            this%rbuffxC(:,:,:,1) = dpFdx
            call this%spectC%fft(this%rbuffxC(:,:,:,1),this%dpF_dxhat)
        end if  

        ! STEP 11: Initialize Wall Model
        if (useWallModelBot) then
            allocate(this%moengWall)
            call this%moengWall%init(this%dz, inputfile, this%gpC, this%sp_gpC, this%rbuffxC, this%rbuffyC, &
                                this%rbuffzC, this%cbuffzC )
            call message(0,"Wall model initialized successfully")
        end if 

        ! STEP 12: Initialize SGS model
        if (this%useSGS) then
            allocate(this%SGSmodel)
            if (allocated(this%moengWall)) then
                call this%sgsModel%init(SGSModelID, this%spectC, this%spectE, this%gpC, this%gpE, this%dx, & 
                    this%dy, this%dz, useDynamicProcedure, useSGSclipping, this%moengWall)
            else
                call this%sgsModel%init(SGSModelID, this%spectC, this%spectE, this%gpC, this%gpE, this%dx, & 
                    this%dy, this%dz, useDynamicProcedure, useSGSclipping)
            end if
            call this%sgsModel%link_pointers(this%nu_SGS, this%c_SGS, this%tauSGS_ij)
            call message(0,"SGS model initialized successfully")
        end if 
        this%max_nuSGS = zero


        ! STEP 13: Set visualization planes for io
        call set_planes_io(this%xplanes, this%yplanes, this%zplanes)


        ! Final Step: Safeguard against unfinished procedures

        call message("IGRID initialized successfully!")
        call message("===========================================================")
    end subroutine


    function getMaxKE(this) result(maxKE)
        class(igrid), intent(inout) :: this
        real(rkind)  :: maxKE
        integer :: i, j, k
        real(rkind) :: keloc

        this%rbuffxC(:,:,:,1) = this%u**2 + this%v**2 + this%wC**2
        maxKE = half*p_maxval(maxval(this%rbuffxC))

    end function

    subroutine interp_What_to_WhatC(this)
        class(igrid), intent(inout), target :: this
        complex(rkind), dimension(:,:,:), pointer :: ybuffC, ybuffE, zbuffC, zbuffE

        ybuffE => this%cbuffyE(:,:,:,1)
        zbuffE => this%cbuffzE(:,:,:,1)
        zbuffC => this%cbuffzC(:,:,:,1)
        ybuffC => this%cbuffyC(:,:,:,1)


        ! Step 1: Transpose what from y -> z
        call transpose_y_to_z(this%what,zbuffE,this%sp_gpE)

        ! Step 2: Interpolate from E -> C
        if (useCompactFD) then
            call this%derW%InterpZ_E2C(zbuffE,zbuffC,this%nxZ, this%nyZ)
        else
            call this%Ops%InterpZ_Edge2Cell(zbuffE,zbuffC)
        end if

        ! Step 3: Transpose back from z -> y
        call transpose_z_to_y(zbuffC,this%whatC,this%sp_gpC)

        ! Step 4: Get wC
        call this%spectC%ifft(this%whatC,this%wC)

        nullify(ybuffE) 
        nullify(zbuffE) 
        nullify(zbuffC) 
        nullify(ybuffC) 


        ! Done !
    end subroutine


    subroutine printDivergence(this)
        class(igrid), intent(inout) :: this

        if (useCompactFD) then
            call this%padepoiss%DivergenceCheck(this%uhat, this%vhat, this%what, this%divergence)
        else
            call this%poiss%DivergenceCheck(this%uhat, this%vhat, this%what, this%divergence)
        end if 
        call message(1, "Domain Maximum Divergence:", p_maxval(this%divergence))
        
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
        call this%filzC%filter3(zbuff3,zbuff4,this%nxZ, this%nyZ)
        call transpose_z_to_y(zbuff4,this%what, this%sp_gpE)

        nullify(zbuff1, zbuff2, zbuff3, zbuff4)
    end subroutine


    !subroutine laplacian(this, f, lapf)
    !    class(igrid),target, intent(inout) :: this
    !    real(rkind), intent(in),  dimension(this%nxp, this%nyp, this%nzp) :: f
    !    real(rkind), intent(out), dimension(this%nxp, this%nyp, this%nzp) :: lapf

    !    lapf = f 
    !end subroutine

    !subroutine gradient(this,f,dfdx,dfdy,dfdz)
    !    class(igrid), intent(inout), target :: this
    !    real(rkind), dimension(this%nxp, this%nyp, this%nzp), intent(in):: f
    !    real(rkind), dimension(this%nxp, this%nyp, this%nzp), intent(out):: dfdx 
    !    real(rkind), dimension(this%nxp, this%nyp, this%nzp), intent(out):: dfdy
    !    real(rkind), dimension(this%nxp, this%nyp, this%nzp), intent(out):: dfdz
    !    
    !    dfdx = f
    !    dfdy = f
    !    dfdz = f
    !end subroutine 

    subroutine destroy(this)
        class(igrid), intent(inout) :: this
        
        nullify(this%u, this%uhat, this%v, this%vhat, this%w, this%what, this%wC)
        deallocate(this%PfieldsC, this%PfieldsE, this%SfieldsC, this%SfieldsE)
        nullify(this%u_rhs, this%v_rhs, this%w_rhs)
        deallocate(this%rhsC, this%rhsE, this%OrhsC, this%OrhsE)
        deallocate(this%duidxj, this%duidxjhat)
        call this%spectC%destroy()
        call this%spectE%destroy()
        deallocate(this%spectC, this%spectE)
        nullify(this%nu_SGS, this%c_SGS, this%tauSGS_ij)
        call this%sgsModel%destroy()
        deallocate(this%sgsModel)
    end subroutine

    subroutine compute_vorticity(this)
        class(igrid), intent(inout), target :: this
        complex(rkind), dimension(:,:,:), pointer :: cbuffz1, cbuffz2
        complex(rkind), dimension(:,:,:), pointer :: cbuffy1
        integer :: i, j, k

        ! Call this subroutine only after all uhat and u fields are the latest
        ! ones

        cbuffz1 => this%cbuffzC(:,:,:,1)
        cbuffz2 => this%cbuffzC(:,:,:,2)
        cbuffy1 => this%cbuffyC(:,:,:,1)

        ! Compute omega_x hat
        call transpose_y_to_z(this%vhat,cbuffz1,this%sp_gpC)
        if (useCompactFD) then
            call this%derV%ddz_C2C(cbuffz1,cbuffz2,this%nxZ,this%nyZ)    
        else
            call this%Ops%ddz_C2C(cbuffz1,cbuffz2,topBC_v,botBC_v)
        end if 
        call transpose_z_to_y(cbuffz2,this%oxhat,this%sp_gpC)
        call this%spectC%mTimes_ik2_oop(this%whatC,cbuffy1) 

        do k = 1,size(this%oxhat,3)
            do j = 1,size(this%oxhat,2)
                do i = 1,size(this%oxhat,1)
                    this%oxhat(i,j,k) = - this%oxhat(i,j,k) + cbuffy1(i,j,k)
                end do 
            end do 
        end do 


        ! Compute omega_y hat
        call transpose_y_to_z(this%uhat,cbuffz1,this%sp_gpC)
        if (useCompactFD) then
            call this%derU%ddz_C2C(cbuffz1,cbuffz2,this%nxZ,this%nyZ)
        else
            call this%Ops%ddz_C2C(cbuffz1,cbuffz2,topBC_u,botBC_u)
        end if 
        call transpose_z_to_y(cbuffz2,this%oyhat,this%sp_gpC)
        call this%spectC%mTimes_ik1_oop(this%whatC,cbuffy1) 
        do k = 1,size(this%oxhat,3)
            do j = 1,size(this%oxhat,2)
                do i = 1,size(this%oxhat,1)
                    this%oyhat(i,j,k) =  this%oyhat(i,j,k) - cbuffy1(i,j,k)
                end do 
            end do 
        end do 

        ! Compute omega_z hat
        !this%ozhat = k1*this%vhat 
        call this%spectC%mTimes_ik1_oop(this%vhat,this%ozhat) 
        call this%spectC%mTimes_ik2_oop(this%uhat,cbuffy1) 
        do k = 1,size(this%oxhat,3)
            do j = 1,size(this%oxhat,2)
                do i = 1,size(this%oxhat,1)
                    this%ozhat(i,j,k) =  this%ozhat(i,j,k) - cbuffy1(i,j,k)
                end do 
            end do 
        end do 
   
        ! Compute ox, oy, oz
        call this%spectC%ifft(this%oxhat,this%ox)
        call this%spectC%ifft(this%oyhat,this%oy)
        call this%spectC%ifft(this%ozhat,this%oz)
        
        nullify(cbuffz1, cbuffz2, cbuffy1)

    end subroutine


    subroutine addNonLinearTerm_Rot(this)
        class(igrid), intent(inout), target :: this
        real(rkind), dimension(:,:,:), pointer :: rtmpx1
        complex(rkind), dimension(:,:,:), pointer :: ctmpz1, ctmpz2
        integer :: i, j, k
        
        ! Assume that the rhs vectors haven't been populated

        rtmpx1 => this%rbuffxC(:,:,:,1)
        ctmpz1 => this%cbuffzC(:,:,:,1)
        ctmpz2 => this%cbuffzE(:,:,:,1)
        
        ! x equation
        do k = 1,size(this%oz,3)
            do j = 1,size(this%oz,2)
                do i = 1,size(this%oz,1)
                    rtmpx1(i,j,k) = this%oz(i,j,k)*this%v(i,j,k) - this%oy(i,j,k)*this%wC(i,j,k)
                end do 
            end do 
        end do 
        call this%spectC%fft(rtmpx1,this%u_rhs)

        ! y equation
        do k = 1,size(this%oz,3)
            do j = 1,size(this%oz,2)
                do i = 1,size(this%oz,1)
                    rtmpx1(i,j,k) = this%ox(i,j,k)*this%wC(i,j,k) - this%oz(i,j,k)*this%u(i,j,k)
                end do 
            end do 
        end do 
        call this%spectC%fft(rtmpx1,this%v_rhs)

        ! z equation 
        do k = 1,size(this%oz,3)
            do j = 1,size(this%oz,2)
                do i = 1,size(this%oz,1)
                    rtmpx1(i,j,k) = this%oy(i,j,k)*this%u(i,j,k) - this%ox(i,j,k)*this%v(i,j,k)
                end do 
            end do 
        end do 
        call this%spectC%fft(rtmpx1,this%wC_rhs)
        
        ! Interpolate w_rhs using wC_rhs
        call transpose_y_to_z(this%wC_rhs,ctmpz1, this%sp_gpC)
        if (useCompactFD) then
            call this%derW%InterpZ_C2E(ctmpz1,ctmpz2,this%nxZ, this%nyZ)
        else
            call this%Ops%InterpZ_Cell2Edge(ctmpz1,ctmpz2,zeroC,zeroC)   
        end if 
        call transpose_z_to_y(ctmpz2,this%w_rhs,this%sp_gpE)

        ! Done 
        nullify(rtmpx1) 
        nullify(ctmpz1) 
        nullify(ctmpz2) 
    end subroutine

    subroutine addNonLinearTerm_SkewSymm(this)
        class(igrid), intent(inout), target :: this
        real(rkind), dimension(:,:,:), pointer :: rtmpx1, rtmpx2
        complex(rkind), dimension(:,:,:), pointer :: ctmpz1, ctmpz2, ctmpz3, ctmpz4
        complex(rkind), dimension(:,:,:), pointer :: ctmpy1, ctmpy2, ctmpy3
        real(rkind), dimension(:,:,:), pointer :: dudx, dudy, dudz
        real(rkind), dimension(:,:,:), pointer :: dvdx, dvdy, dvdz
        real(rkind), dimension(:,:,:), pointer :: dwdx, dwdy, dwdz

        dudx => this%duidxj(:,:,:,1); dudy => this%duidxj(:,:,:,2); dudz => this%duidxj(:,:,:,3); 
        dvdx => this%duidxj(:,:,:,4); dvdy => this%duidxj(:,:,:,5); dvdz => this%duidxj(:,:,:,6); 
        dwdx => this%duidxj(:,:,:,7); dwdy => this%duidxj(:,:,:,8); dwdz => this%duidxj(:,:,:,9); 


        rtmpx1 => this%rbuffxC(:,:,:,1); rtmpx2 => this%rbuffxE(:,:,:,1)
        
        ctmpz1 => this%cbuffzC(:,:,:,1); ctmpz2 => this%cbuffzC(:,:,:,2)
        ctmpz3 => this%cbuffzE(:,:,:,1); ctmpz4 => this%cbuffzE(:,:,:,2)
        
        ctmpy1 => this%cbuffyC(:,:,:,1); ctmpy2 => this%cbuffyE(:,:,:,1)
        ctmpy3 => this%cbuffyC(:,:,:,2)

        ! Advection Form Terms
        rtmpx1 = -this%u*dudx
        call this%spectC%fft(rtmpx1,this%u_rhs)
        rtmpx1 = -this%v*dudy
        call this%spectC%fft(rtmpx1,ctmpy1)
        this%u_rhs = this%u_rhs + ctmpy1
        rtmpx1 = -this%wC*dudz
        call this%spectC%fft(rtmpx1,ctmpy1)
        this%u_rhs = this%u_rhs + ctmpy1

        rtmpx1 = -this%u*dvdx
        call this%spectC%fft(rtmpx1,this%v_rhs)
        rtmpx1 = -this%v*dvdy
        call this%spectC%fft(rtmpx1,ctmpy1)
        this%v_rhs = this%v_rhs + ctmpy1
        rtmpx1 = -this%wC*dvdz
        call this%spectC%fft(rtmpx1,ctmpy1)
        this%v_rhs = this%v_rhs + ctmpy1
        
        rtmpx1 = -this%u*dwdx
        call this%spectC%fft(rtmpx1,ctmpy1)
        call transpose_y_to_z(ctmpy1,ctmpz1,this%sp_gpC)
        if (useCompactFD) then
            call this%derW%interpZ_C2E(ctmpz1,ctmpz3,size(ctmpz1,1),size(ctmpz1,2))
        else
            call this%Ops%interpZ_Cell2Edge(ctmpz1,ctmpz3,zeroC,zeroC)
        end if 
        call transpose_z_to_y(ctmpz3,this%w_rhs,this%sp_gpE)

        rtmpx1 = -this%v*dwdy
        call this%spectC%fft(rtmpx1,ctmpy1)
        call transpose_y_to_z(ctmpy1,ctmpz1,this%sp_gpC)
        if (useCompactFD) then
            call this%derW%interpZ_C2E(ctmpz1,ctmpz3,size(ctmpz1,1),size(ctmpz1,2))
        else
            call this%Ops%interpZ_Cell2Edge(ctmpz1,ctmpz3,zeroC,zeroC)
        end if 
        call transpose_z_to_y(ctmpz3,ctmpy2,this%sp_gpE)
        this%w_rhs = this%w_rhs + ctmpy2

        rtmpx1 = -this%wC*dwdz
        call this%spectC%fft(rtmpx1,ctmpy1)
        call transpose_y_to_z(ctmpy1,ctmpz1,this%sp_gpC)
        if (useCompactFD) then
            call this%derW%interpZ_C2E(ctmpz1,ctmpz3,size(ctmpz1,1),size(ctmpz1,2))
        else
            call this%Ops%interpZ_Cell2Edge(ctmpz1,ctmpz3,zeroC,zeroC)
        end if 
        call transpose_z_to_y(ctmpz3,ctmpy2,this%sp_gpE)
        this%w_rhs = this%w_rhs + ctmpy2


        this%u_rhs = half*this%u_rhs
        this%v_rhs = half*this%v_rhs
        this%w_rhs = half*this%w_rhs

        ! Conservative Form Terms        
        rtmpx1 = -this%u*this%u
        call this%spectC%fft(rtmpx1,ctmpy1)
        call this%spectC%mtimes_ik1_ip(ctmpy1)
        this%u_rhs = this%u_rhs + half*ctmpy1
        
        rtmpx1 = -this%u*this%v
        call this%spectC%fft(rtmpx1,ctmpy1)
        call this%spectC%mtimes_ik2_oop(ctmpy1,ctmpy3)
        this%u_rhs = this%u_rhs + half*ctmpy3
        call this%spectC%mtimes_ik1_oop(ctmpy1,ctmpy3)
        this%v_rhs = this%v_rhs + half*ctmpy3

        rtmpx1 = -this%u*this%wC
        call this%spectC%fft(rtmpx1,ctmpy1)
        call transpose_y_to_z(ctmpy1,ctmpz1,this%sp_gpC)
        
        if (useCompactFD) then
            call this%derW%ddz_C2C(ctmpz1,ctmpz2,size(ctmpz1,1),size(ctmpz1,2))
        else
            call this%Ops%ddz_C2C(ctmpz1,ctmpz2,topBC_w,botBC_w)
        end if 
        
        call transpose_z_to_y(ctmpz2,ctmpy1,this%sp_gpC)
        this%u_rhs = this%u_rhs + half*ctmpy1
        
        if (useCompactFD) then
            call this%derW%interpZ_C2E(ctmpz1,ctmpz3,size(ctmpz1,1),size(ctmpz1,2))
        else
            call this%Ops%interpZ_Cell2Edge(ctmpz1,ctmpz3,zeroC,zeroC)
        end if 
        
        call transpose_z_to_y(ctmpz3,ctmpy2,this%sp_gpE)
        call this%spectE%mtimes_ik1_ip(ctmpy2)
        this%w_rhs = this%w_rhs + half*ctmpy2

        rtmpx1 = -this%v*this%v
        call this%spectC%fft(rtmpx1,ctmpy1)
        call this%spectC%mtimes_ik2_ip(ctmpy1)
        this%v_rhs = this%v_rhs + half*ctmpy1


        rtmpx1 = -this%v*this%wC
        call this%spectC%fft(rtmpx1,ctmpy1)
        call transpose_y_to_z(ctmpy1,ctmpz1,this%sp_gpC)
        
        if (useCompactFD) then
            call this%derW%ddz_C2C(ctmpz1,ctmpz2,size(ctmpz1,1),size(ctmpz1,2))
        else
            call this%Ops%ddz_C2C(ctmpz1,ctmpz2,topBC_w,botBC_w)
        end if 
        
        call transpose_z_to_y(ctmpz2,ctmpy1,this%sp_gpC)
        this%v_rhs = this%v_rhs + half*ctmpy1

        if (useCompactFD) then
            call this%derW%InterpZ_C2E(ctmpz1,ctmpz3,size(ctmpz1,1),size(ctmpz1,2))
        else
            call this%Ops%interpZ_Cell2Edge(ctmpz1,ctmpz3,zeroC,zeroC)
        end if 

        call transpose_z_to_y(ctmpz3,ctmpy2,this%sp_gpE)
        call this%spectE%mtimes_ik2_ip(ctmpy2)
        this%w_rhs = this%w_rhs + half*ctmpy2

        rtmpx1 = -this%wC*this%wC
        call this%spectC%fft(rtmpx1,ctmpy1)
        call transpose_y_to_z(ctmpy1,ctmpz1,this%sp_gpC)
        
        if (useCompactFD) then
            call this%derWW%ddz_C2E(ctmpz1,ctmpz3,size(ctmpz1,1),size(ctmpz1,2))
        else
            call this%Ops%ddz_C2E(ctmpz1,ctmpz3,.true.,.true.)
        end if 

        call transpose_z_to_y(ctmpz3,ctmpy2,this%sp_gpE)
        this%w_rhs = this%w_rhs + half*ctmpy2


        nullify( dudx, dudy, dudz) 
        nullify( dvdx, dvdy, dvdz)
        nullify( dwdx, dwdy, dwdz)
    end subroutine


    subroutine addCoriolisTerm(this)
        class(igrid), intent(inout) :: this
        integer :: i, j, k

        ! u equation 
        do k = 1,size(this%u_rhs,3)
            do j = 1,size(this%u_rhs,2)
                do i = 1,size(this%u_rhs,1)
                    this%u_rhs(i,j,k) = this%u_rhs(i,j,k) + this%vhat(i,j,k)/this%Ro
                end do 
            end do 
        end do 
        

        ! v equation 
        do k = 1,size(this%u_rhs,3)
            do j = 1,size(this%u_rhs,2)
                do i = 1,size(this%u_rhs,1)
                    this%v_rhs(i,j,k) = this%v_rhs(i,j,k) +  (this%GxHat(i,j,k) - this%uhat(i,j,k))/this%Ro
                end do 
            end do 
        end do 

        ! w equation 
        ! Do nothing 
    end subroutine  


    subroutine addExtraForcingTerm(this)
        class(igrid), intent(inout) :: this
       
        this%u_rhs = this%u_rhs + this%dpF_dxhat

    end subroutine

    subroutine addViscousTerm(this)
        class(igrid), intent(inout), target :: this
        complex(rkind), dimension(:,:,:), pointer :: cytmp1, cztmp1, cztmp2
        complex(rkind), dimension(:,:,:), pointer :: cztmp3, cztmp4, cytmp2
        integer :: i, j, k
        
        cytmp1 => this%cbuffyC(:,:,:,1)
        cytmp2 => this%cbuffyE(:,:,:,1)
        cztmp1 => this%cbuffzC(:,:,:,1)
        cztmp2 => this%cbuffzC(:,:,:,2)
        cztmp3 => this%cbuffzE(:,:,:,1)
        cztmp4 => this%cbuffzE(:,:,:,2)
        
        ! u equation 
        call transpose_y_to_z(this%uhat, cztmp1,this%sp_gpC)
        if (useCompactFD) then
            call this%derU%d2dz2_C2C(cztmp1,cztmp2,this%nxZ,this%nyZ)
        else
            call this%Ops%d2dz2_C2C(cztmp1,cztmp2,topBC_u,botBC_u)
        end if 
        
        call transpose_z_to_y(cztmp2,cytmp1,this%sp_gpC)
        do k = 1,size(this%u_rhs,3)
            do j = 1,size(this%u_rhs,2)
                do i = 1,size(this%u_rhs,1)
                    this%u_rhs(i,j,k) = this%u_rhs(i,j,k) + (one/this%Re)*cytmp1(i,j,k) &
                             - (one/this%Re)*this%spectC%kabs_sq(i,j,k) *this%uhat(i,j,k)
                end do 
            end do 
        end do 

        ! v equation 
        call transpose_y_to_z(this%vhat, cztmp1,this%sp_gpC)
        if (useCompactFD) then
            call this%derV%d2dz2_C2C(cztmp1,cztmp2,this%nxZ,this%nyZ)
        else
            call this%Ops%d2dz2_C2C(cztmp1,cztmp2,topBC_v,botBC_v)
        end if 
        call transpose_z_to_y(cztmp2,cytmp1,this%sp_gpC)
        do k = 1,size(this%v_rhs,3)
            do j = 1,size(this%v_rhs,2)
                do i = 1,size(this%v_rhs,1)
                    this%v_rhs(i,j,k) = this%v_rhs(i,j,k) + (one/this%Re)*cytmp1(i,j,k) &
                             - (one/this%Re)*this%spectC%kabs_sq(i,j,k) *this%vhat(i,j,k)
                end do 
            end do 
        end do 
        
        ! w equation
        call transpose_y_to_z(this%what, cztmp3,this%sp_gpE)
        if (useCompactFD) then
            call this%derW%d2dz2_C2C(cztmp3,cztmp4,this%nxZ,this%nyZ)
        else
            call this%Ops%d2dz2_E2E(cztmp3,cztmp4,topBC_w,botBC_w)
        end if 

        call transpose_z_to_y(cztmp4,cytmp2,this%sp_gpC)
        do k = 1,size(this%w_rhs,3)
            do j = 1,size(this%w_rhs,2)
                do i = 1,size(this%w_rhs,1)
                    this%w_rhs(i,j,k) = this%w_rhs(i,j,k) + (one/this%Re)*cytmp2(i,j,k) &
                             - (one/this%Re)*this%spectE%kabs_sq(i,j,k) *this%what(i,j,k)
                end do 
            end do 
        end do 

        nullify(cytmp1,cytmp2,cztmp1,cztmp2,cztmp3,cztmp4)
    end subroutine 

    subroutine AdamsBashforth(this)
        class(igrid), intent(inout) :: this
        integer :: i, j, k

        ! Step 1: Non Linear Term 
        select case(AdvectionForm)
        case (1) ! Rotational Form
            call this%AddNonLinearTerm_Rot()
        case (2) ! SkewSymm Form
            call this%AddNonLinearTerm_SkewSymm()
        end select  

        ! Step 2: Coriolis Term
        if (this%useCoriolis) then
            call this%AddCoriolisTerm()
        end if 
       
        if (this%useExtraForcing) then
            call this%addExtraForcingTerm()
        end if 

        ! Step 3a: Viscous Term
        if (.not. this%isInviscid) then
            call this%AddViscousTerm()
        end if 

        ! Step 3b: SGS Viscous Term
        if (this%useSGS) then
            call this%SGSmodel%getRHS_SGS(this%duidxj, this%duidxjhat, this%u_rhs, &
                                          this%v_rhs , this%w_rhs    , this%uhat , &
                                          this%vhat  , this%whatC    , this%u    , &
                                          this%v     , this%wC       , this%max_nuSGS)

            !! IMPORTANT: duidxj, u, v and wC are all corrupted if SGS was initialized to use the
            !! Dynamic Procedure. DON'T USE duidxj again within this time step.
            !! Make the SGS call at the very end, just before the time
            !! advancement.
        end if 

        ! Step 4: Time Step 
        if (this%step == 0) then
            do k = 1,size(this%uhat,3)
                do j = 1,size(this%uhat,2)
                    do i = 1,size(this%uhat,1)
                        this%uhat(i,j,k) = this%uhat(i,j,k) + this%dt*this%u_rhs(i,j,k)
                    end do 
                end do 
            end do  
            do k = 1,size(this%vhat,3)
                do j = 1,size(this%vhat,2)
                    do i = 1,size(this%vhat,1)
                        this%vhat(i,j,k) = this%vhat(i,j,k) + this%dt*this%v_rhs(i,j,k)
                    end do 
                end do 
            end do  
            do k = 1,size(this%what,3)
                do j = 1,size(this%what,2)
                    do i = 1,size(this%what,1)
                        this%what(i,j,k) = this%what(i,j,k) + this%dt*this%w_rhs(i,j,k)
                    end do 
                end do 
            end do  
        else
            do k = 1,size(this%uhat,3)
                do j = 1,size(this%uhat,2)
                    do i = 1,size(this%uhat,1)
                        this%uhat(i,j,k) = this%uhat(i,j,k) + this%dtby2*(three*this%u_rhs(i,j,k) - &
                                                                this%u_Orhs(i,j,k))
                    end do 
                end do 
            end do  
            do k = 1,size(this%vhat,3)
                do j = 1,size(this%vhat,2)
                    do i = 1,size(this%vhat,1)
                        this%vhat(i,j,k) = this%vhat(i,j,k) + this%dtby2*(three*this%v_rhs(i,j,k) - &
                                                                this%v_Orhs(i,j,k))
                    end do 
                end do 
            end do  
            do k = 1,size(this%what,3)
                do j = 1,size(this%what,2)
                    do i = 1,size(this%what,1)
                        this%what(i,j,k) = this%what(i,j,k) + this%dtby2*(three*this%w_rhs(i,j,k) - &
                                                                this%w_Orhs(i,j,k))
                    end do 
                end do 
            end do 
        end if 
        

        ! Step 5: Dealias 
        call this%spectC%dealias(this%uhat)
        call this%spectC%dealias(this%vhat)
        call this%spectE%dealias(this%what)
        if (this%useVerticalFilter) then
            call this%ApplyCompactFilter()
        end if 

        
        ! Step 6: Pressure projection
        if (useCompactFD) then
            call this%padepoiss%PressureProjection(this%uhat,this%vhat,this%what)
            call this%padepoiss%DivergenceCheck(this%uhat, this%vhat, this%what, this%divergence, .true.) 
        else
            call this%poiss%PressureProjNP(this%uhat,this%vhat,this%what)
            call this%poiss%DivergenceCheck(this%uhat, this%vhat, this%what, this%divergence) 
        end if 

        ! Step 7: Take it back to physical fields
        call this%spectC%ifft(this%uhat,this%u)
        call this%spectC%ifft(this%vhat,this%v)
        call this%spectE%ifft(this%what,this%w)

        ! STEP 8: Interpolate the cell center values of w
        call this%interp_wHat_to_wHatC()


        ! STEP 9: Compute vorticity and duidxj 
        if (AdvectionForm == 1) then
            call this%compute_Vorticity()
        end if 
        if ((AdvectionForm == 2) .or. (this%useSGS)) then
            call this%compute_duidxj()
        end if 
        
        ! STEP 10: Copy the RHS for using during next time step 
        do k = 1,size(this%uhat,3)
            do j = 1,size(this%uhat,2)
                do i = 1,size(this%uhat,1)
                    this%u_Orhs(i,j,k) = this%u_rhs(i,j,k)
                end do 
            end do
        end do  
        do k = 1,size(this%vhat,3)
            do j = 1,size(this%vhat,2)
                do i = 1,size(this%vhat,1)
                    this%v_Orhs(i,j,k) = this%v_rhs(i,j,k)
                end do 
            end do
        end do  
        do k = 1,size(this%what,3)
            do j = 1,size(this%what,2)
                do i = 1,size(this%what,1)
                    this%w_Orhs(i,j,k) = this%w_rhs(i,j,k)
                end do 
            end do
        end do  


    end subroutine

    subroutine compute_duidxj(this)
        class(igrid), intent(inout), target :: this
        complex(rkind), dimension(:,:,:), pointer :: ctmpz1, ctmpz2
        complex(rkind), dimension(:,:,:), pointer :: ctmpz3, ctmpz4
        complex(rkind), dimension(:,:,:), pointer :: ctmpy1, ctmpy2
        real(rkind),    dimension(:,:,:), pointer :: dudx, dudy, dudz
        real(rkind),    dimension(:,:,:), pointer :: dvdx, dvdy, dvdz
        real(rkind),    dimension(:,:,:), pointer :: dwdx, dwdy, dwdz
        complex(rkind), dimension(:,:,:), pointer :: dudxH, dudyH, dudzH 
        complex(rkind), dimension(:,:,:), pointer :: dvdxH, dvdyH, dvdzH
        complex(rkind), dimension(:,:,:), pointer :: dwdxH, dwdyH, dwdzH

        dudx => this%duidxj(:,:,:,1); dudy => this%duidxj(:,:,:,2); dudz => this%duidxj(:,:,:,3); 
        dvdx => this%duidxj(:,:,:,4); dvdy => this%duidxj(:,:,:,5); dvdz => this%duidxj(:,:,:,6); 
        dwdx => this%duidxj(:,:,:,7); dwdy => this%duidxj(:,:,:,8); dwdz => this%duidxj(:,:,:,9); 

        dudxH => this%duidxjhat(:,:,:,1); dudyH => this%duidxjhat(:,:,:,2); dudzH => this%duidxjhat(:,:,:,3); 
        dvdxH => this%duidxjhat(:,:,:,4); dvdyH => this%duidxjhat(:,:,:,5); dvdzH => this%duidxjhat(:,:,:,6); 
        dwdxH => this%duidxjhat(:,:,:,7); dwdyH => this%duidxjhat(:,:,:,8); dwdzH => this%duidxjhat(:,:,:,9); 
        
        ctmpz1 => this%cbuffzC(:,:,:,1); ctmpz2 => this%cbuffzE(:,:,:,1); 
        ctmpz3 => this%cbuffzC(:,:,:,2); ctmpz4 => this%cbuffzE(:,:,:,2)
        
        ctmpy1 => this%cbuffyC(:,:,:,1); ctmpy2 => this%cbuffyE(:,:,:,1)

        call this%spectC%mTimes_ik1_oop(this%uhat,dudxH)
        call this%spectC%ifft(dudxH,dudx)

        call this%spectC%mTimes_ik2_oop(this%uhat,dudyH)
        call this%spectC%ifft(dudyH,dudy)

        call this%spectC%mTimes_ik1_oop(this%vhat,dvdxH)
        call this%spectC%ifft(dvdxH,dvdx)

        call this%spectC%mTimes_ik2_oop(this%vhat,dvdyH)
        call this%spectC%ifft(dvdyH,dvdy)
        
        call this%spectC%mTimes_ik1_oop(this%whatC,dwdxH)
        call this%spectC%ifft(dwdxH,dwdx)

        call this%spectC%mTimes_ik2_oop(this%whatC,dwdyH)
        call this%spectC%ifft(dwdyH,dwdy)

        call transpose_y_to_z(this%uhat,ctmpz1, this%sp_gpC)
        if (useCompactFD) then
            call this%derU%ddz_C2C(ctmpz1,ctmpz3,size(ctmpz1,1),size(ctmpz1,2))
        else
            call this%Ops%ddz_C2C(ctmpz1,ctmpz3,topBC_u,botBC_u)
        end if 
        call transpose_z_to_y(ctmpz3,dudzH,this%sp_gpC)
        call this%spectC%ifft(dudzH,dudz)

        call transpose_y_to_z(this%vhat,ctmpz1, this%sp_gpC)
        if (useCompactFD) then
            call this%derV%ddz_C2C(ctmpz1,ctmpz3,size(ctmpz1,1),size(ctmpz1,2))
        else    
            call this%Ops%ddz_C2C(ctmpz1,ctmpz3,topBC_v,botBC_v) 
        end if 
        call transpose_z_to_y(ctmpz3,dvdzH,this%sp_gpC)
        call this%spectC%ifft(dvdzH,dvdz)

        call transpose_y_to_z(this%what,ctmpz2, this%sp_gpE)
        if (useCompactFD) then
            call this%derW%ddz_E2C(ctmpz2,ctmpz3,size(ctmpz2,1),size(ctmpz2,2))
        else
            call this%Ops%ddz_E2C(ctmpz2,ctmpz3) 
        end if 
        call transpose_z_to_y(ctmpz3,dwdzH,this%sp_gpC)
        call this%spectC%ifft(dwdzH,dwdz)
        
        nullify( dudx, dudy, dudz) 
        nullify( dvdx, dvdy, dvdz)
        nullify( dwdx, dwdy, dwdz)
        nullify( dudxH, dudyH, dudzH) 
        nullify( dvdxH, dvdyH, dvdzH)
        nullify( dwdxH, dwdyH, dwdzH)
        nullify( ctmpy1,ctmpy2,ctmpz1, ctmpz2, ctmpz3,ctmpz4)

    end subroutine


    subroutine readRestartFile(this, tid, rid)
        use decomp_2d_io
        use mpi
        use exits, only: message
        class(igrid), intent(inout) :: this
        integer, intent(in) :: tid, rid
        character(len=clen) :: tempname, fname
        integer :: ierr

        write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",rid, "_u.",tid
        fname = this%InputDir(:len_trim(this%InputDir))//"/"//trim(tempname)
        call decomp_2d_read_one(1,this%u,fname)

        write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",rid, "_v.",tid
        fname = this%InputDir(:len_trim(this%InputDir))//"/"//trim(tempname)
        call decomp_2d_read_one(1,this%v,fname)

        write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",rid, "_w.",tid
        fname = this%InputDir(:len_trim(this%InputDir))//"/"//trim(tempname)
        call decomp_2d_read_one(1,this%w,fname)

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
        class(igrid), intent(in) :: this
        character(len=clen) :: tempname, fname
        integer :: ierr

        write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",this%runID, "_u.",this%step
        fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
        call decomp_2d_write_one(1,this%u,fname)

        write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",this%runID, "_v.",this%step
        fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
        call decomp_2d_write_one(1,this%v,fname)

        write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",this%runID, "_w.",this%step
        fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
        call decomp_2d_write_one(1,this%w,fname)

        write(tempname,"(A7,A4,I2.2,A6,I6.6)") "RESTART", "_Run",this%runID, "_info.",this%step
        fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)

        OPEN(UNIT=10, FILE=trim(fname))
        write(10,"(100g15.5)") this%tsim
        close(10)

        call mpi_barrier(mpi_comm_world, ierr)
        call message(1, "Just Dumped a RESTART file")

    end subroutine 

    !! STATISTICS !!

    subroutine init_stats( this)
        class(igrid), intent(inout), target :: this
        type(decomp_info), pointer  :: gpC

        gpC => this%gpC
        this%tidSUM = 0

        allocate(this%zStats2dump(this%nz,24))
        allocate(this%runningSum(this%nz,24))
        allocate(this%TemporalMnNOW(this%nz,24))

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

        this%runningSum = zero
        nullify(gpC)
    end subroutine

    subroutine compute_stats(this)
        class(igrid), intent(inout), target :: this
        type(decomp_info), pointer :: gpC
        real(rkind), dimension(:,:,:), pointer :: rbuff1, rbuff2, rbuff3, rbuff4, rbuff5, rbuff6

        rbuff1 => this%rbuffxC(:,:,:,1); rbuff2 => this%rbuffyC(:,:,:,1);
        rbuff3 => this%rbuffzC(:,:,:,1); 
        rbuff4 => this%rbuffzC(:,:,:,2); 
        rbuff5 => this%rbuffzC(:,:,:,3); 
        rbuff6 => this%rbuffzC(:,:,:,4); 
        !rbuff7 => this%rbuffzC(:,:,:,5); 
        gpC => this%gpC

        this%tidSUM = this%tidSUM + 1

        ! Compute u - mean 
        call transpose_x_to_y(this%u,rbuff2,this%gpC)
        call transpose_y_to_z(rbuff2,rbuff3,this%gpC)
        call this%compute_z_mean(rbuff3, this%u_mean)

        ! Compute v - mean 
        call transpose_x_to_y(this%v,rbuff2,this%gpC)
        call transpose_y_to_z(rbuff2,rbuff4,this%gpC)
        call this%compute_z_mean(rbuff4, this%v_mean)

        ! Compute w - mean 
        call transpose_x_to_y(this%wC,rbuff2,this%gpC)
        call transpose_y_to_z(rbuff2,rbuff5,this%gpC)
        call this%compute_z_mean(rbuff5, this%w_mean)

        ! uu mean
        rbuff6 = rbuff3*rbuff3
        call this%compute_z_mean(rbuff6, this%uu_mean)

        ! uv mean
        rbuff6 = rbuff3*rbuff4
        call this%compute_z_mean(rbuff6, this%uv_mean)

        ! uw mean
        rbuff6 = rbuff3*rbuff5
        call this%compute_z_mean(rbuff6, this%uw_mean)

        ! vv mean 
        rbuff6 = rbuff4*rbuff4
        call this%compute_z_mean(rbuff6, this%vv_mean)

        ! vw mean 
        rbuff6 = rbuff4*rbuff5
        call this%compute_z_mean(rbuff6, this%vw_mean)

        ! ww mean 
        rbuff6 = rbuff5*rbuff5
        call this%compute_z_mean(rbuff6, this%ww_mean)

        if (this%useSGS) then
            ! tau_11
            call transpose_x_to_y(this%tauSGS_ij(:,:,:,1),rbuff2,this%gpC)
            call transpose_y_to_z(rbuff2,rbuff3,this%gpC)
            call this%compute_z_mean(rbuff3, this%tau11_mean)

            ! tau_12
            call transpose_x_to_y(this%tauSGS_ij(:,:,:,2),rbuff2,this%gpC)
            call transpose_y_to_z(rbuff2,rbuff3,this%gpC)
            call this%compute_z_mean(rbuff3, this%tau12_mean)

            ! tau_13
            call transpose_x_to_y(this%tauSGS_ij(:,:,:,3),rbuff2,this%gpC)
            call transpose_y_to_z(rbuff2,rbuff3,this%gpC)
            call this%compute_z_mean(rbuff3, this%tau13_mean)

            ! tau_22
            call transpose_x_to_y(this%tauSGS_ij(:,:,:,4),rbuff2,this%gpC)
            call transpose_y_to_z(rbuff2,rbuff3,this%gpC)
            call this%compute_z_mean(rbuff3, this%tau22_mean)

            ! tau_23
            call transpose_x_to_y(this%tauSGS_ij(:,:,:,5),rbuff2,this%gpC)
            call transpose_y_to_z(rbuff2,rbuff3,this%gpC)
            call this%compute_z_mean(rbuff3, this%tau23_mean)

            ! tau_33
            call transpose_x_to_y(this%tauSGS_ij(:,:,:,6),rbuff2,this%gpC)
            call transpose_y_to_z(rbuff2,rbuff3,this%gpC)
            call this%compute_z_mean(rbuff3, this%tau33_mean)


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
        class(igrid), intent(inout), target :: this
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
        class(igrid), intent(in), target :: this
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
        class(igrid), intent(in), target :: this
        real(rkind), dimension(:,:,:), intent(in) :: arr_in
        real(rkind), dimension(:), intent(out) :: vec_out
        integer :: k

        do k = 1,size(arr_in,3)
            vec_out(k) = P_SUM(sum(arr_in(:,:,k)))/(real(this%nx,rkind)*real(this%ny,rkind))
        end do 

    end subroutine

    subroutine finalize_stats(this)
        class(igrid), intent(inout) :: this
        nullify(this%u_mean, this%v_mean, this%w_mean, this%uu_mean, this%uv_mean, this%uw_mean, this%vv_mean, this%vw_mean, this%ww_mean)
        nullify(this%tau11_mean, this%tau12_mean, this%tau13_mean, this%tau22_mean, this%tau23_mean, this%tau33_mean)
        nullify(this%S11_mean, this%S12_mean, this%S13_mean, this%S22_mean, this%S23_mean, this%S33_mean)
        nullify(this%sgsdissp, this%viscdissp, this%sgscoeff_mean)
        deallocate(this%zStats2dump, this%runningSum, this%TemporalMnNOW)
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






end module 
