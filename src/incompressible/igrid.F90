module IncompressibleGridNP
    use kind_parameters, only: rkind, clen
    use constants, only: imi, zero,one,two,three,half 
    use GridMod, only: grid
    use gridtools, only: alloc_buffs, destroy_buffs
    use igrid_hooks, only: meshgen, initfields_stagg, getforcing
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
    use sigmaSGSmod, only: sigmasgs

    implicit none

    private
    public :: igrid 

    logical, parameter :: useCompactFD = .true. 
    integer, parameter :: AdvectionForm = 2 

    integer, parameter :: no_slip = 1, slip = 2
    complex(rkind), parameter :: zeroC = zero + imi*zero 


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
        type(sigmaSGS), allocatable :: SGS

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

        complex(rkind), dimension(:,:,:), allocatable :: dPf_dxhat

        real(rkind) :: max_nuSGS

        contains
            procedure :: init
            procedure :: destroy
            procedure :: printDivergence 
            procedure :: laplacian
            procedure :: gradient
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
        logical :: useSGS = .false. 
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
        integer :: restartFile_TID
        integer :: restartFile_RID
        logical :: useRestartFile = .false. 
        logical :: useWallModelTop = .false., useWallModelBot = .false.
        logical :: isInviscid = .false., useVerticalFilter = .true.  


        namelist /INPUT/       nx, ny, nz, tstop, dt, CFL, nsteps, &
                                              inputdir, outputdir, &
                                  periodicx, periodicy, periodicz, &
                                                       prow, pcol, &
                                        t_restartDump, t_dataDump
        namelist /IINPUT/  Re, useSGS, runID, Pr, useCoriolis, & 
                                tid_statsDump, useExtraForcing, &
                                time_startDumping, topWall, botWall, &
                                useRestartFile, restartFile_TID, restartFile_RID, &
                                useWallModelTop, useWallModelBot, isInviscid, &
                                useVerticalFilter

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

        this%useVerticalFilter = useVerticalFilter

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
        else
            call message("WARNING: No Top BCs provided. Using defaults found in igrid.F90")
        end if 
        
        if (botWall == slip) then
            botBC_u = .true.; botBC_v = .true.
            call message(1, "BotWall BC set to: SLIP")
        elseif (botWall == no_slip) then
            botBC_u = .false.; botBC_v = .false.
            call message(1, "BotWall BC set to: NO_SLIP")
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
            call this%derU%init (nz  , this%dz, topBC_u, botBC_u, useWallModelTop, useWallModelBot) 
            call this%derV%init (nz  , this%dz, topBC_v, botBC_v, useWallModelTop, useWallModelBot)  
            call this%derW%init (nz  , this%dz, topBC_w, botBC_w)
            call this%derWW%init(nz , this%dz, .true., .true.)
        else
            allocate(this%Ops)
            call this%Ops%init(this%gpC,this%gpE,0,this%dx,this%dy,this%dz,this%spectC%spectdecomp,this%spectE%spectdecomp)
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
        allocate(this%rbuffzC(this%gpC%zsz(1),this%gpC%zsz(2),this%gpC%zsz(3),2))

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
            call this%spectC%alloc_r2c_out(this%GxHat)
            this%rbuffxC(:,:,:,1) = this%Gx
            call this%spectC%fft(this%rbuffxC(:,:,:,1),this%Gxhat)
        end if

        ! STEP 10b: Compute additional forcing (channel)
        if (this%useExtraForcing) then
            call getForcing(dpFdx)
            call message(0," Turning on aditional forcing")
            call message(1," dP_dx = ", dpFdx)
            call this%spectC%alloc_r2c_out(this%dpF_dxhat)
            this%rbuffxC(:,:,:,1) = dpFdx
            call this%spectC%fft(this%rbuffxC(:,:,:,1),this%dpF_dxhat)
        end if  

        ! STEP 12: Initialize SGS model
        if (this%useSGS) then
            allocate(this%SGS)
            call this%sgs%init(this%spectC, this%spectE, this%gpC, this%gpE, this%dx, this%dy, this%dz)
            call message(0,"SGS model initialized successfully")
        end if 
        this%max_nuSGS = zero

        ! Final Step: Safeguard against unfinished procedures
        if ((.not.useCompactFD) .and. (AdvectionForm == 2)) then
            call GracefulExit("Skew Symmetric Form is not allowed with 2nd order & 
                & FD. Use 6th order Compact FD scheme instead",213)
        end if 

        call message("IGRID initialized successfully!")
        call message("===========================================================")
    end subroutine


    pure function getMaxKE(this) result(maxKE)
        class(igrid), intent(in) :: this
        real(rkind)  :: maxKE
        integer :: i, j, k
        real(rkind) :: keloc

        maxKE = zero
        do k = 1,size(this%u,3)
            do j = 1,size(this%u,2)
                do i = 1,size(this%u,1)
                    keloc = this%u(i,j,k)**2 
                    keloc = keloc + this%v(i,j,k)**2
                    keloc = keloc + this%w(i,j,k)**2
                    maxKE = max(maxKE,keloc)
                end do 
            end do 
        end do 

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


    subroutine laplacian(this, f, lapf)
        class(igrid),target, intent(inout) :: this
        real(rkind), intent(in),  dimension(this%nxp, this%nyp, this%nzp) :: f
        real(rkind), intent(out), dimension(this%nxp, this%nyp, this%nzp) :: lapf

        lapf = f 
    end subroutine

    subroutine gradient(this,f,dfdx,dfdy,dfdz)
        class(igrid), intent(inout), target :: this
        real(rkind), dimension(this%nxp, this%nyp, this%nzp), intent(in):: f
        real(rkind), dimension(this%nxp, this%nyp, this%nzp), intent(out):: dfdx 
        real(rkind), dimension(this%nxp, this%nyp, this%nzp), intent(out):: dfdy
        real(rkind), dimension(this%nxp, this%nyp, this%nzp), intent(out):: dfdz
        
        dfdx = f
        dfdy = f
        dfdz = f
    end subroutine 

    subroutine destroy(this)
        class(igrid), intent(inout) :: this
        
        nullify(this%u, this%uhat, this%v, this%vhat, this%w, this%what, this%wC)
        deallocate(this%PfieldsC, this%PfieldsE, this%SfieldsC, this%SfieldsE)
        nullify(this%u_rhs, this%v_rhs, this%w_rhs)
        deallocate(this%rhsC, this%rhsE, this%OrhsC, this%OrhsE)
        call this%spectC%destroy()
        call this%spectE%destroy()
        deallocate(this%spectC, this%spectE)
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
        call this%derW%interpZ_C2E(ctmpz1,ctmpz3,size(ctmpz1,1),size(ctmpz1,2))
        call transpose_z_to_y(ctmpz3,this%w_rhs,this%sp_gpE)

        rtmpx1 = -this%v*dwdy
        call this%spectC%fft(rtmpx1,ctmpy1)
        call transpose_y_to_z(ctmpy1,ctmpz1,this%sp_gpC)
        call this%derW%interpZ_C2E(ctmpz1,ctmpz3,size(ctmpz1,1),size(ctmpz1,2))
        call transpose_z_to_y(ctmpz3,ctmpy2,this%sp_gpE)
        this%w_rhs = this%w_rhs + ctmpy2

        rtmpx1 = -this%wC*dwdz
        call this%spectC%fft(rtmpx1,ctmpy1)
        call transpose_y_to_z(ctmpy1,ctmpz1,this%sp_gpC)
        call this%derW%interpZ_C2E(ctmpz1,ctmpz3,size(ctmpz1,1),size(ctmpz1,2))
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
        call this%derW%ddz_C2C(ctmpz1,ctmpz2,size(ctmpz1,1),size(ctmpz1,2))
        call transpose_z_to_y(ctmpz2,ctmpy1,this%sp_gpC)
        this%u_rhs = this%u_rhs + half*ctmpy1
        call this%derW%InterpZ_C2E(ctmpz1,ctmpz3,size(ctmpz1,1),size(ctmpz1,2))
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
        call this%derW%ddz_C2C(ctmpz1,ctmpz2,size(ctmpz1,1),size(ctmpz1,2))
        call transpose_z_to_y(ctmpz2,ctmpy1,this%sp_gpC)
        this%v_rhs = this%v_rhs + half*ctmpy1
        call this%derW%InterpZ_C2E(ctmpz1,ctmpz3,size(ctmpz1,1),size(ctmpz1,2))
        call transpose_z_to_y(ctmpz3,ctmpy2,this%sp_gpE)
        call this%spectE%mtimes_ik2_ip(ctmpy2)
        this%w_rhs = this%w_rhs + half*ctmpy2

        rtmpx1 = -this%wC*this%wC
        call this%spectC%fft(rtmpx1,ctmpy1)
        call transpose_y_to_z(ctmpy1,ctmpz1,this%sp_gpC)
        call this%derWW%ddz_C2E(ctmpz1,ctmpz3,size(ctmpz1,1),size(ctmpz1,2))
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
            call this%SGS%getRHS_SGS(this%duidxj,this%u_rhs,this%v_rhs,this%w_rhs, this%max_nuSGS)
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

        dudx => this%duidxj(:,:,:,1); dudy => this%duidxj(:,:,:,2); dudz => this%duidxj(:,:,:,3); 
        dvdx => this%duidxj(:,:,:,4); dvdy => this%duidxj(:,:,:,5); dvdz => this%duidxj(:,:,:,6); 
        dwdx => this%duidxj(:,:,:,7); dwdy => this%duidxj(:,:,:,8); dwdz => this%duidxj(:,:,:,9); 

        ctmpz1 => this%cbuffzC(:,:,:,1)
        ctmpz2 => this%cbuffzE(:,:,:,1)
        ctmpz3 => this%cbuffzC(:,:,:,2)
        ctmpz4 => this%cbuffzE(:,:,:,2)
        
        ctmpy1 => this%cbuffyC(:,:,:,1)
        ctmpy2 => this%cbuffyE(:,:,:,1)


        call this%spectC%mTimes_ik1_oop(this%uhat,ctmpy1)
        call this%spectC%ifft(ctmpy1,dudx)

        call this%spectC%mTimes_ik2_oop(this%uhat,ctmpy1)
        call this%spectC%ifft(ctmpy1,dudy)

        call this%spectC%mTimes_ik1_oop(this%vhat,ctmpy1)
        call this%spectC%ifft(ctmpy1,dvdx)

        call this%spectC%mTimes_ik2_oop(this%vhat,ctmpy1)
        call this%spectC%ifft(ctmpy1,dvdy)
        
        call this%spectC%mTimes_ik1_oop(this%whatC,ctmpy1)
        call this%spectC%ifft(ctmpy1,dwdx)

        call this%spectC%mTimes_ik2_oop(this%whatC,ctmpy1)
        call this%spectC%ifft(ctmpy1,dwdy)

        call transpose_y_to_z(this%uhat,ctmpz1, this%sp_gpC)
        if (useCompactFD) then
            call this%derU%ddz_C2C(ctmpz1,ctmpz3,size(ctmpz1,1),size(ctmpz1,2))
        else
            call this%Ops%ddz_C2C(ctmpz1,ctmpz3,topBC_u,botBC_u) 
        end if 
        call transpose_z_to_y(ctmpz3,ctmpy1,this%sp_gpC)
        call this%spectC%ifft(ctmpy1,dudz)

        call transpose_y_to_z(this%vhat,ctmpz1, this%sp_gpC)
        if (useCompactFD) then
            call this%derV%ddz_C2C(ctmpz1,ctmpz3,size(ctmpz1,1),size(ctmpz1,2))
        else    
            call this%Ops%ddz_C2C(ctmpz1,ctmpz3,topBC_v,botBC_v) 
        end if 
        call transpose_z_to_y(ctmpz3,ctmpy1,this%sp_gpC)
        call this%spectC%ifft(ctmpy1,dvdz)

        call transpose_y_to_z(this%what,ctmpz2, this%sp_gpE)
        if (useCompactFD) then
            call this%derW%ddz_E2C(ctmpz2,ctmpz3,size(ctmpz2,1),size(ctmpz2,2))
        else
            call this%Ops%ddz_E2C(ctmpz2,ctmpz3) 
        end if 
        call transpose_z_to_y(ctmpz3,ctmpy1,this%sp_gpC)
        call this%spectC%ifft(ctmpy1,dwdz)
        
        nullify( dudx, dudy, dudz) 
        nullify( dvdx, dvdy, dvdz)
        nullify( dwdx, dwdy, dwdz)
        nullify(ctmpy1,ctmpy2,ctmpz1, ctmpz2, ctmpz3,ctmpz4)

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

end module 

