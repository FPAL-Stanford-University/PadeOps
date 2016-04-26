module IncompressibleGridWallM
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
    use sgsmod, only: sgs
    use wallmodelMod, only: wallmodel
    use numerics, only: use3by2rule 
    implicit none

    private
    public :: igridWallM 

    real(rkind), parameter :: kappa = 0.41d0
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
        type(staggOps), allocatable :: Ops
        type(sgs), allocatable :: sgsmodel
        type(wallmodel), allocatable :: moengWall

        real(rkind), dimension(:,:,:,:), allocatable :: PfieldsC
        real(rkind), dimension(:,:,:,:), allocatable :: PfieldsE

        complex(rkind), dimension(:,:,:,:), allocatable :: SfieldsC
        complex(rkind), dimension(:,:,:,:), allocatable :: SfieldsE


        type(poisson), allocatable :: poiss
        type(padepoisson), allocatable :: padepoiss
        real(rkind), dimension(:,:,:), allocatable :: divergence

        real(rkind), dimension(:,:,:), pointer :: u, v, wC, w, uE, vE
        real(rkind), dimension(:,:,:), pointer :: ox,oy,oz
        complex(rkind), dimension(:,:,:), pointer :: uhat, vhat, whatC, what


        real(rkind), dimension(:,:,:,:), allocatable, public :: rbuffxC, rbuffyC, rbuffzC
        real(rkind), dimension(:,:,:,:), allocatable :: rbuffxE, rbuffyE, rbuffzE
        
        complex(rkind), dimension(:,:,:,:), allocatable :: cbuffyC, cbuffzC
        complex(rkind), dimension(:,:,:,:), allocatable :: cbuffyE, cbuffzE

        complex(rkind), dimension(:,:,:,:), allocatable :: rhsC, rhsE, OrhsC, OrhsE 
        real(rkind), dimension(:,:,:,:), allocatable :: duidxjC, duidxjE 
        complex(rkind), dimension(:,:,:,:), allocatable :: duidxjChat
        complex(rkind), dimension(:,:,:), pointer:: u_rhs, v_rhs, wC_rhs, w_rhs 
        complex(rkind), dimension(:,:,:), pointer:: u_Orhs, v_Orhs, w_Orhs
            
        real(rkind) :: Re, Gx, Gy, Gz, dtby2, meanfact
        complex(rkind), dimension(:,:,:), allocatable :: GxHat 
        real(rkind) :: Ro = 1.d5

        integer :: nxZ, nyZ
        integer :: tid_statsDump
        real(rkind) :: time_startDumping 
        
        integer :: runID
        logical :: useCoriolis = .true. 
        logical :: useExtraForcing = .false.
        logical :: useSGS = .false. 
        logical :: useVerticalFilter = .false. 
        logical :: useDynamicProcedure 
        logical :: useCFL = .false.  
        logical :: dumpPlanes = .false.

        complex(rkind), dimension(:,:,:), allocatable :: dPf_dxhat

        real(rkind) :: max_nuSGS

        real(rkind) :: z0, ustar, Umn, Vmn, Uspmn, dtOld, dtRat
        real(rkind), dimension(:,:,:), allocatable :: filteredSpeedSq
        integer :: wallMType 

        ! Statistics to compute 
        real(rkind), dimension(:,:), allocatable :: zStats2dump, runningSum, TemporalMnNOW
        real(rkind), dimension(:), pointer :: u_mean, v_mean, w_mean, uu_mean, uv_mean, uw_mean, vv_mean, vw_mean, ww_mean
        real(rkind), dimension(:), pointer :: tau11_mean, tau12_mean, tau13_mean, tau22_mean, tau23_mean, tau33_mean
        real(rkind), dimension(:), pointer :: S11_mean, S12_mean, S13_mean, S22_mean, S23_mean, S33_mean
        real(rkind), dimension(:), pointer :: viscdissp, sgsdissp, sgscoeff_mean, PhiM
        integer :: tidSUM


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
            procedure, private :: addNonLinearTerm_Rot
            procedure, private :: addCoriolisTerm
            procedure, private :: addExtraForcingTerm 
            procedure, private :: compute_and_bcast_surface_Mn
            procedure          :: dumpRestartFile
            procedure, private :: readRestartFile
            procedure, private :: compute_z_mean 
            procedure, private :: compute_z_fluct
            procedure, private :: compute_deltaT
            procedure, private :: getfilteredSpeedSqAtWall
            procedure          :: dump_stats
            procedure          :: compute_stats 
            procedure          :: finalize_stats
            procedure          :: dump_planes
            procedure          :: dumpFullField 
    end type

contains 

    subroutine init(this,inputfile)
        class(igridWallM), intent(inout), target :: this        
        character(len=clen), intent(in) :: inputfile 
        integer :: nx, ny, nz
        character(len=clen) :: outputdir
        character(len=clen) :: inputdir
        integer :: prow = 0, pcol = 0 
        logical :: useSGS = .true., useDynamicProcedure = .false.
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
        real(rkind) :: Pr = 0.7_rkind 
        integer :: topWall = slip
        logical :: useCoriolis = .true. 
        logical :: useExtraForcing = .false.
        real(rkind) :: dpFdx = zero
        integer :: restartFile_TID = 1
        integer :: restartFile_RID = 1
        logical :: useRestartFile = .false. 
        logical :: isInviscid = .false., useVerticalFilter = .true.  
        integer :: SGSModelID = 1, WallMType = 0
        real(rkind) :: z0 = 1.d-4
        logical :: dumpPlanes = .false. 

        namelist /INPUT/       nx, ny, nz, tstop, dt, CFL, nsteps, &
                                              inputdir, outputdir, &
                                                       prow, pcol, &
                                        t_restartDump, t_dataDump
        namelist /IINPUT/  Re, useDynamicProcedure,runID, Pr, useCoriolis, & 
                                tid_statsDump, useExtraForcing, useSGSclipping, &
                                time_startDumping, topWall, &
                                useRestartFile, restartFile_TID, restartFile_RID, &
                                isInviscid, dumpPlanes, &
                                useVerticalFilter, SGSModelID

        namelist /WALLMODEL/ z0, wallMType 

        ! STEP 1: READ INPUT 
        ioUnit = 11
        open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
        read(unit=ioUnit, NML=INPUT)
        read(unit=ioUnit, NML=IINPUT)
        read(unit=ioUnit, NML=WALLMODEL)
        close(ioUnit)
      
        this%nx = nx; this%ny = ny; this%nz = nz
        this%meanfact = one/(real(nx,rkind)*real(ny,rkind))
        this%dt = dt; this%dtby2 = dt/two ; this%z0 = z0 ; this%Re = Re
        this%outputdir = outputdir; this%inputdir = inputdir 
        this%WallMtype = WallMType
        this%runID = runID; this%tstop = tstop; this%t_dataDump = t_dataDump
        this%CFL = CFL; this%dumpPlanes = dumpPlanes
        if (this%CFL > zero) this%useCFL = .true. 
        if ((this%CFL < zero) .and. (this%dt < zero)) then
            call GracefulExit("Both CFL and dt cannot be negative. Have you &
            & specified either one of these in the input file?", 124)
        end if 

        this%t_restartDump = t_restartDump; this%tid_statsDump = tid_statsDump
        this%time_startDumping = time_startDumping ; this%useCoriolis = useCoriolis 
        this%useExtraForcing = useExtraForcing; this%useSGS = useSGS 
        this%useVerticalFilter = useVerticalFilter
        this%useDynamicProcedure = useDynamicProcedure

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
            call message("WARNING: No Top BCs provided. Using defaults found in igridWallM.F90")
        end if 
        
        

        ! STEP 3: GENERATE MESH (CELL CENTERED) 
        if ( allocated(this%mesh) ) deallocate(this%mesh) 
        allocate(this%mesh(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3),3))
        call meshgen(this%gpC, this%dx, this%dy, &
            this%dz, this%mesh) ! <-- this procedure is part of user defined HOOKS
        
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
        allocate(this%Ops)
        call this%Ops%init(this%gpC,this%gpE,0,this%dx,this%dy,this%dz,this%spectC%spectdecomp, &
                    this%spectE%spectdecomp, .false., .false.)
        
        ! STEP 6: ALLOCATE MEMORY FOR FIELD ARRAYS
        allocate(this%PfieldsC(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3),6))
        allocate(this%divergence(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)))
        allocate(this%duidxjC(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3),9))
        allocate(this%duidxjE(this%gpE%xsz(1),this%gpE%xsz(2),this%gpE%xsz(3),4))
        allocate(this%PfieldsE(this%gpE%xsz(1),this%gpE%xsz(2),this%gpE%xsz(3),3))
        call this%spectC%alloc_r2c_out(this%SfieldsC,3); call this%spectC%alloc_r2c_out(this%duidxjChat,9)
        call this%spectC%alloc_r2c_out(this%rhsC,3); call this%spectC%alloc_r2c_out(this%OrhsC,2)
        call this%spectE%alloc_r2c_out(this%rhsE,1); call this%spectE%alloc_r2c_out(this%OrhsE,1)
        call this%spectE%alloc_r2c_out(this%SfieldsE,1)
        
        this%u => this%PfieldsC(:,:,:,1) ; this%v => this%PfieldsC(:,:,:,2) ; this%wC => this%PfieldsC(:,:,:,3) 
        this%w => this%PfieldsE(:,:,:,1) ; this%uE => this%PfieldsE(:,:,:,2) ; this%vE => this%PfieldsE(:,:,:,3) 
        
        this%uhat => this%SfieldsC(:,:,:,1); this%vhat => this%SfieldsC(:,:,:,2); 
        this%whatC => this%SfieldsC(:,:,:,3); this%what => this%SfieldsE(:,:,:,1)

        this%ox => this%PfieldsC(:,:,:,4); this%oy => this%PfieldsC(:,:,:,5); this%oz => this%PfieldsC(:,:,:,6)

        this%u_rhs => this%rhsC(:,:,:,1); this%v_rhs => this%rhsC(:,:,:,2); this%w_rhs => this%rhsE(:,:,:,1)

        this%u_Orhs => this%OrhsC(:,:,:,1); this%v_Orhs => this%OrhsC(:,:,:,2); this%w_Orhs => this%OrhsE(:,:,:,1)

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
        allocate(this%filteredSpeedSq(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)))
        this%nxZ = size(this%cbuffzE,1); this%nyZ = size(this%cbuffzE,2)

        ! STEP 6: ALLOCATE/INITIALIZE THE POISSON DERIVED TYPE
        allocate(this%poiss)
        call this%poiss%init(this%spectC,.false.,this%dx,this%dy,this%dz,this%Ops,this%spectE, .true.)  
               
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
        
        this%Gx = u_g; this%Gy = zero; this%Gz = zero

        call this%spectC%fft(this%u,this%uhat)   
        call this%spectC%fft(this%v,this%vhat)   
        call this%spectE%fft(this%w,this%what)   
      
        ! Dealias and filter before projection
        call this%spectC%dealias(this%uhat)
        call this%spectC%dealias(this%vhat)
        call this%spectE%dealias(this%what)

        ! Pressure projection
        call this%poiss%PressureProjNP(this%uhat,this%vhat,this%what)
        call this%poiss%DivergenceCheck(this%uhat, this%vhat, this%what, this%divergence)

        ! Take it back to physical fields
        call this%spectC%ifft(this%uhat,this%u)
        call this%spectC%ifft(this%vhat,this%v)
        call this%spectE%ifft(this%what,this%w)
       
        ! STEP 8: Interpolate the cell center values of w
        call this%compute_and_bcast_surface_Mn()
        call this%interp_PrimitiveVars()
        call message(1,"Max KE:",P_MAXVAL(this%getMaxKE()))
     
        ! STEP 9: Compute duidxj
         call this%compute_duidxj()
       

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

        ! STEP 11: Initialize SGS model
        if (this%useSGS) then
            allocate(this%SGSmodel)
            call this%SGSmodel%init(SGSModelID, this%spectC, this%spectE, this%gpC, this%gpE, this%dx, & 
                this%dy, this%dz, useDynamicProcedure, useSGSclipping, this%mesh(:,:,:,3),1, this%z0, .true., WallMType, .true.)
            call this%sgsModel%link_pointers(this%nu_SGS, this%c_SGS, this%tauSGS_ij, this%tau13, this%tau23)
            call message(0,"SGS model initialized successfully")
        end if 
        this%max_nuSGS = zero


        ! STEP 12: Set visualization planes for io
        call set_planes_io(this%xplanes, this%yplanes, this%zplanes)


        ! STEP 13: Compute the timestep
        call this%compute_deltaT()
        this%dtOld = this%dt
        this%dtRat = one 

        ! Final Step: Safeguard against unfinished procedures

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
        call this%Ops%InterpZ_Edge2Cell(zbuffE,zbuffC)
        call transpose_z_to_y(zbuffC,this%whatC,this%sp_gpC)
        call this%spectC%ifft(this%whatC,this%wC)

        ! Step 2: Interpolate u -> uE
        call transpose_y_to_z(this%uhat,zbuffC,this%sp_gpC)
        call this%Ops%InterpZ_Cell2Edge(zbuffC,zbuffE,zeroC,zeroC)
        
        zbuffE(:,:,1) = (three/two)*zbuffC(:,:,1) - half*zbuffC(:,:,2)
        zbuffE(:,:,this%nz + 1) = zbuffC(:,:,this%nz)

        call transpose_z_to_y(zbuffE,ybuffE,this%sp_gpE)
        call this%spectE%ifft(ybuffE,this%uE)
        
        ! Step 3: Interpolate v -> vE
        call transpose_y_to_z(this%vhat,zbuffC,this%sp_gpC)
        call this%Ops%InterpZ_Cell2Edge(zbuffC,zbuffE,zeroC,zeroC)
        
        zbuffE(:,:,1) = (three/two)*zbuffC(:,:,1) - half*zbuffC(:,:,2)
        zbuffE(:,:,this%nz + 1) = zbuffC(:,:,this%nz)
        
        call transpose_z_to_y(zbuffE,ybuffE,this%sp_gpE)
        call this%spectE%ifft(ybuffE,this%vE)

    end subroutine


    subroutine printDivergence(this)
        class(igridWallM), intent(inout) :: this

        call this%poiss%DivergenceCheck(this%uhat, this%vhat, this%what, this%divergence)
        call message(1, "Domain Maximum Divergence:", p_maxval(this%divergence))
        
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
        
        dudx  => this%duidxjC(:,:,:,1); dudy  => this%duidxjC(:,:,:,2); dudzC => this%duidxjC(:,:,:,3); 
        dvdx  => this%duidxjC(:,:,:,4); dvdy  => this%duidxjC(:,:,:,5); dvdzC => this%duidxjC(:,:,:,6); 
        dwdxC => this%duidxjC(:,:,:,7); dwdyC => this%duidxjC(:,:,:,8); dwdz  => this%duidxjC(:,:,:,9); 

        dwdx => this%duidxjE(:,:,:,1); dwdy => this%duidxjE(:,:,:,2);
        dudz => this%duidxjE(:,:,:,3); dvdz => this%duidxjE(:,:,:,4);

        T1C => this%rbuffxC(:,:,:,1); T2C => this%rbuffxC(:,:,:,2)
        T1E => this%rbuffxE(:,:,:,1); T2E => this%rbuffxE(:,:,:,2)
       
        ! Step 1: u - equation rhs
        T1C = dvdx - dudy
        T1C = T1c*this%v
        T2C = dwdxC - dudzC
        T2C = T2C*this%wC
        T1C = T1C + T2C
        call this%spectC%fft(T1C,this%u_rhs)

        ! Step 2: v - equation rhs
        T1C = dudy - dvdx
        T1C = T1C*this%u
        T2C = dwdyC - dvdzC
        T2C = T2C*this%wC
        T1C = T1C + T2C
        call this%spectC%fft(T1C,this%v_rhs)

        ! Step 3: w - equation
        T1E = dudz - dwdx
        T1E = T1E*this%uE
        T2E = dvdz - dwdy
        T2E = T2E*this%vE
        T1E = T1E + T2E
        call this%spectE%fft(T1E,this%w_rhs)

    end subroutine


    subroutine addCoriolisTerm(this)
        class(igridWallM), intent(inout) :: this
        ! u equation 
        this%u_rhs = this%u_rhs + this%vhat/this%Ro
        ! v equation 
        this%v_rhs = this%v_rhs +  (this%GxHat - this%uhat)/this%Ro
    end subroutine  


    subroutine addExtraForcingTerm(this)
        class(igridWallM), intent(inout) :: this
        !if (this%spectC%carryingZeroK) then
        !    this%dpF_dxhat(1,1,:) = cmplx(this%ustar*this%ustar*this%nx*this%ny,zero)
        !end if
        this%u_rhs = this%u_rhs + this%dpF_dxhat
    end subroutine

    
    subroutine AdamsBashforth(this)
        class(igridWallM), intent(inout) :: this
        real(rkind) :: abf1, abf2

        ! Step 0: Compute TimeStep 
        call this%compute_deltaT
        this%dtRat = this%dt/this%dtOld

        ! Step 1: Non Linear Term 
        call this%AddNonLinearTerm_Rot()

        ! Step 2: Coriolis Term
        if (this%useCoriolis) then
            call this%AddCoriolisTerm()
        end if 
       
        if (this%useExtraForcing) then
            call this%addExtraForcingTerm()
        end if 

        ! Step 3b: SGS Viscous Term
        if (this%useSGS) then
            call this%SGSmodel%getRHS_SGS_WallM(this%duidxjC, this%duidxjE  , this%duidxjChat,& 
                                                this%u_rhs  , this%v_rhs    , this%w_rhs     ,&
                                                this%uhat   , this%vhat     , this%whatC     ,&
                                                this%u      , this%v        , this%wC        ,&
                                                this%ustar  , this%Umn      , this%Vmn       ,&
                                                this%Uspmn  , this%filteredSpeedSq, this%max_nuSGS)

            !! IMPORTANT: duidxjC, u, v and wC are all corrupted if SGS was initialized to use the
            !! Dynamic Procedure. DON'T USE duidxjC again within this time step.
            !! Make the SGS call at the very end, just before the time
            !! advancement.
        end if 
        

        ! Step 4: Time Step 
        if (this%step == 0) then
            this%uhat = this%uhat + this%dt*this%u_rhs 
            this%vhat = this%vhat + this%dt*this%v_rhs 
            this%what = this%what + this%dt*this%w_rhs 
        else
            abf1 = (one + half*this%dtRat)*this%dt
            abf2 = -half*this%dtRat*this%dt
            this%uhat = this%uhat + abf1*this%u_rhs + abf2*this%u_Orhs
            this%vhat = this%vhat + abf1*this%v_rhs + abf2*this%v_Orhs
            this%what = this%what + abf1*this%w_rhs + abf2*this%w_Orhs
        end if 
        

        ! Step 5: Dealias
        call this%spectC%dealias(this%uhat)
        call this%spectC%dealias(this%vhat)
        call this%spectE%dealias(this%what)
        
        ! Step 6: Pressure projection
        call this%poiss%PressureProjNP(this%uhat,this%vhat,this%what)
        call this%poiss%DivergenceCheck(this%uhat, this%vhat, this%what, this%divergence) 

        ! Step 7: Take it back to physical fields
        call this%spectC%ifft(this%uhat,this%u)
        call this%spectC%ifft(this%vhat,this%v)
        call this%spectE%ifft(this%what,this%w)

        ! STEP 8: Interpolate the cell center values of w
        call this%compute_and_bcast_surface_Mn()
        call this%interp_PrimitiveVars()

        ! STEP 9: Compute duidxjC 
        call this%compute_duidxj()
        
        ! STEP 10: Copy the RHS for using during next time step 
        this%u_Orhs = this%u_rhs
        this%v_Orhs = this%v_rhs
        this%w_Orhs = this%w_rhs
        this%dtOld = this%dt

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
        call this%Ops%ddz_E2C(ctmpz2,ctmpz1)
        call transpose_z_to_y(ctmpz1,dwdzH,this%sp_gpC)
        call this%spectC%ifft(dwdzH,dwdz)


        ! Compute dudz, dvdz 
        call transpose_y_to_z(this%uhat,ctmpz1,this%sp_gpC)
        call this%Ops%ddz_C2E(ctmpz1,ctmpz2,topBC_u,.true.)
            

        !dudz_dzby2 = ctmpz1(:,:,1)/((this%dz/two)*log(this%dz/two/this%z0))
        !call this%spectE%SurfaceFilter_ip(dudz_dzby2)
        !dudz_dzby2 = (this%Umn/this%Uspmn) * dudz_dzby2

        !ctmpz2(:,:,1) = two*dudz_dzby2 - ctmpz2(:,:,2)
        ctmpz2(:,:,1) = (two*ctmpz2(:,:,2) - ctmpz2(:,:,3))
        call transpose_z_to_y(ctmpz2,ctmpy2,this%sp_gpE)
        call this%spectE%ifft(ctmpy2,dudz)
        call this%Ops%InterpZ_Edge2Cell(ctmpz2,ctmpz1)
        call transpose_z_to_y(ctmpz1,dudzH,this%sp_gpC)
        call this%spectC%ifft(dudzH,dudzC)

        call transpose_y_to_z(this%vhat,ctmpz1,this%sp_gpC)
        call this%Ops%ddz_C2E(ctmpz1,ctmpz2,topBC_v,.true.)

        !dvdz_dzby2 = dudz_dzby2 * this%Vmn/this%Umn
        !ctmpz2(:,:,1) = two*dvdz_dzby2 - ctmpz2(:,:,2)
        ctmpz2(:,:,1) = (two*ctmpz2(:,:,2) - ctmpz2(:,:,3))
        call transpose_z_to_y(ctmpz2,ctmpy2,this%sp_gpE)
        call this%spectE%ifft(ctmpy2,dvdz)
        call this%Ops%InterpZ_Edge2Cell(ctmpz2,ctmpz1)
        call transpose_z_to_y(ctmpz1,dvdzH,this%sp_gpC)
        call this%spectC%ifft(dvdzH,dvdzC)

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
        end if
        call mpi_bcast(this%Umn,1,mpirkind,0,mpi_comm_world,ierr)
        call mpi_bcast(this%Vmn,1,mpirkind,0,mpi_comm_world,ierr)
        call mpi_bcast(this%Uspmn,1,mpirkind,0,mpi_comm_world,ierr)

        call this%getfilteredSpeedSqAtWall()
        ! Compute USTAR and Umn at z = dz/2
        this%ustar = this%Umn*kappa/log(this%dz/two/this%z0)
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
        real(rkind), dimension(:,:,:), pointer :: rbuff1, rbuff2, rbuff3E, rbuff2E, rbuff3, rbuff4, rbuff5, rbuff6

        rbuff1  => this%rbuffxC(:,:,:,1); rbuff2  => this%rbuffyC(:,:,:,1);
        rbuff2E => this%rbuffyE(:,:,:,1); rbuff3E => this%rbuffzE(:,:,:,1);
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
            call transpose_x_to_y(this%tau13,rbuff2E,this%gpE)
            call transpose_y_to_z(rbuff2E,rbuff3E,this%gpE)
            rbuff3E(:,:,1) = -(this%ustar**2)
            call this%Ops%InterpZ_Edge2Cell(rbuff3E,rbuff3)
            call this%compute_z_mean(rbuff3, this%tau13_mean)

            ! tau_22
            call transpose_x_to_y(this%tauSGS_ij(:,:,:,4),rbuff2,this%gpC)
            call transpose_y_to_z(rbuff2,rbuff3,this%gpC)
            call this%compute_z_mean(rbuff3, this%tau22_mean)

            ! tau_23
            call transpose_x_to_y(this%tau23,rbuff2E,this%gpE)
            call transpose_y_to_z(rbuff2E,rbuff3E,this%gpE)
            rbuff3E(:,:,1) = -(this%ustar**2)*this%Vmn/this%Umn
            call this%Ops%InterpZ_Edge2Cell(rbuff3E,rbuff3)
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
