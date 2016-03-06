module IncompressibleGridNP
    use kind_parameters, only: rkind, clen
    use constants, only: imi, zero,one,two,three,half 
    use GridMod, only: grid
    use gridtools, only: alloc_buffs, destroy_buffs
    use hooks, only: meshgen, initfields_stagg
    use decomp_2d
    use StaggOpsMod, only: staggOps  
    use exits, only: GracefulExit, message
    use spectralMod, only: spectral  
    use PoissonMod, only: poisson
    use mpi 
    use reductions, only: p_maxval

    implicit none

    private
    public :: igrid 

    complex(rkind), parameter :: zeroC = zero + imi*zero 

    integer, parameter :: no_slip = 1, slip = 2

    ! Allow non-zero value (isEven) 
    logical :: topBC_u = .true.  , topBC_v = .true. , topBC_w = .false.
    logical :: botBC_u = .false. , botBC_v = .false., botBC_w = .false. 

    integer, parameter :: AdvectionForm = 3 

    type, extends(grid) :: igrid
        
        type(decomp_info), allocatable :: gpC, gpE
        type(decomp_info), pointer :: Sp_gpC, Sp_gpE
        type(spectral), allocatable :: spectE, spectC
        type(staggOps), allocatable :: Ops

        real(rkind), dimension(:,:,:,:), allocatable :: PfieldsC
        real(rkind), dimension(:,:,:,:), allocatable :: PfieldsE

        complex(rkind), dimension(:,:,:,:), allocatable :: SfieldsC
        complex(rkind), dimension(:,:,:,:), allocatable :: SfieldsE


        type(poisson), allocatable :: poiss
        real(rkind), dimension(:,:,:), allocatable :: divergence

        real(rkind), dimension(:,:,:), pointer :: u, v, wC, w
        real(rkind), dimension(:,:,:), pointer :: ox,oy,oz
        complex(rkind), dimension(:,:,:), pointer :: uhat, vhat, whatC, what
        complex(rkind), dimension(:,:,:), pointer :: oxhat, oyhat, ozhat

        real(rkind), dimension(:,:,:,:), allocatable, public :: rbuffxC, rbuffyC, rbuffzC
        real(rkind), dimension(:,:,:,:), allocatable :: rbuffxE, rbuffyE, rbuffzE
        
        complex(rkind), dimension(:,:,:,:), allocatable :: cbuffyC, cbuffzC
        complex(rkind), dimension(:,:,:,:), allocatable :: cbuffyE, cbuffzE

        complex(rkind), dimension(:,:,:,:), allocatable :: duidxj, rhsC, rhsE, OrhsC, OrhsE 
        complex(rkind), dimension(:,:,:), pointer:: u_rhs, v_rhs, wC_rhs, w_rhs 
        complex(rkind), dimension(:,:,:), pointer:: u_Orhs, v_Orhs, w_Orhs
            
        real(rkind) :: nu0, Gx, Gy, Gz, fCor, dtby2

        integer :: tid_statsDump
        real(rkind) :: time_startDumping 
        
        integer :: runID
        logical :: useCoriolis = .true. 
        contains
            procedure :: init
            procedure :: destroy
            procedure :: printDivergence 
            procedure :: laplacian
            procedure :: gradient
            procedure :: AdamsBashforth
            procedure, private :: interp_wHat_to_wHatC
            procedure, private :: interp_w_to_wC
            procedure :: compute_vorticity
            procedure, private :: compute_duidxj
            procedure, private :: addNonLinearTerm_Rot
            procedure, private :: addNonLinearTerm_Cnsrv
            procedure, private :: addNonLinearTerm_SkewSymm 
            procedure, private :: addCoriolisTerm
            procedure, private :: addViscousTerm 
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
        real(rkind) :: nu = 0.02_rkind
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
        namelist /INPUT/       nx, ny, nz, tstop, dt, CFL, nsteps, &
                                              inputdir, outputdir, &
                                  periodicx, periodicy, periodicz, &
                                                       prow, pcol, &
                                        t_restartDump, t_dataDump
        namelist /IINPUT/  nu, useSGS, runID, Pr, useCoriolis, & 
                                tid_statsDump, &
                                time_startDumping, topWall, botWall 

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

        this%outputdir = outputdir 
        
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
        this%step = 0
        this%tstop = tstop 
        this%t_dataDump = t_dataDump

        this%tid_statsDump = tid_statsDump
        this%time_startDumping = time_startDumping 
        this%useCoriolis = useCoriolis 

        select case(topWall)
        case(slip)
            topBC_u = .false.
            topBC_v = .false.
        case(no_slip)
            topBC_u = .true.
            topBC_v = .true.
        case default 
            call GracefulExit("Incorrect choice for topWall. Only two choices &
            & allowed: 1 (no slip) or 2 (slip)",101) 
        end select 

        select case(botWall)
        case(slip)
            botBC_u = .false.
            botBC_v = .false.
        case(no_slip)
            botBC_u = .true.
            botBC_v = .true.
        case default 
            call GracefulExit("Incorrect choice for botWall. Only two choices & 
            & allowed: 1 (no slip) or 2 (slip)",101) 
        end select 

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
        allocate(this%Ops)
        call this%Ops%init(this%gpC,this%gpE,0,this%dx,this%dy,this%dz,this%spectC%spectdecomp,this%spectE%spectdecomp)
        
        ! STEP 6: ALLOCATE MEMORY FOR FIELD ARRAYS
        allocate(this%PfieldsC(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3),6))
        allocate(this%divergence(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)))
        allocate(this%PfieldsE(this%gpE%xsz(1),this%gpE%xsz(2),this%gpE%xsz(3),1))
        call this%spectC%alloc_r2c_out(this%SfieldsC,6)
        call this%spectC%alloc_r2c_out(this%duidxj,9)
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

        ! STEP 6: ALLOCATE/INITIALIZE THE POISSON DERIVED TYPE 
        allocate(this%poiss)
        call this%poiss%init(this%spectC,.false.,this%dx,this%dy,this%dz,this%Ops,this%spectE)  
        
        ! STEP 7: INITIALIZE THE FIELDS 
        call initfields_stagg(this%gpC, this%gpE, this%dx, this%dy, this%dz, &
            inputfile, this%mesh, this%PfieldsC, this%PfieldsE, u_g, this%fcor)! <-- this procedure is part of user defined HOOKS

        this%nu0 = nu
        this%Gx = u_g
        this%Gy = zero
        this%Gz = zero

        call this%spectC%fft(this%u,this%uhat)   
        call this%spectC%fft(this%v,this%vhat)   
        call this%spectE%fft(this%w,this%what)   

        ! Dealias before projection
        this%uhat = this%uhat*this%spectC%Gdealias 
        this%vhat = this%vhat*this%spectC%Gdealias 
        this%what = this%what*this%spectE%Gdealias 


        ! Pressure projection
        call this%poiss%PressureProjNP(this%uhat,this%vhat,this%what)

        ! Take it back to physical fields
        call this%spectC%ifft(this%uhat,this%u)
        call this%spectC%ifft(this%vhat,this%v)
        call this%spectE%ifft(this%what,this%w)


        ! STEP 8: ALLOCATE STORAGE FOR BUFFERS AND DUIDXJHAT
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


        ! STEP 9: Interpolate the cell center values of w
        call this%interp_w_to_wC()
        call this%interp_wHat_to_wHatC()
        call message(1,"Max KE:",P_MAXVAL(half*(this%u**2 + this%v**2 + this%wC**2)))


        ! STEP 10: Compute Vorticity 
        call this%compute_Vorticity()

    end subroutine

    subroutine interp_W_to_WC(this)
        class(igrid), intent(inout), target :: this
        real(rkind), dimension(:,:,:), pointer :: ybuffC, ybuffE, zbuffC, zbuffE

        ybuffE => this%rbuffyE(:,:,:,1)
        zbuffE => this%rbuffzE(:,:,:,1)
        zbuffC => this%rbuffzC(:,:,:,1)
        ybuffC => this%rbuffyC(:,:,:,1)


        ! Step 1: Transpose w from x -> z
        call transpose_x_to_y(this%w,ybuffE,this%gpE)
        call transpose_y_to_z(ybuffE,zbuffE,this%gpE)

        ! Step 2: Interpolate from E -> C
        call this%Ops%InterpZ_Edge2Cell(zbuffE,zbuffC)

        ! Step 3: Transpose back from z -> x
        call transpose_z_to_y(zbuffC,ybuffC,this%gpC)
        call transpose_y_to_x(ybuffC,this%wC,this%gpC)

        nullify(ybuffE) 
        nullify(zbuffE) 
        nullify(zbuffC) 
        nullify(ybuffC) 
        ! Done !
    end subroutine

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
        call this%Ops%InterpZ_Edge2Cell(zbuffE,zbuffC)

        ! Step 3: Transpose back from z -> y
        call transpose_z_to_y(zbuffC,this%whatC,this%sp_gpC)

        nullify(ybuffE) 
        nullify(zbuffE) 
        nullify(zbuffC) 
        nullify(ybuffC) 


        ! Done !
    end subroutine


    subroutine printDivergence(this)
        class(igrid), intent(inout) :: this
        call this%poiss%DivergenceCheck(this%uhat, this%vhat, this%what, this%divergence)
        call message(1, "Domain Maximum Divergence:", p_maxval(this%divergence))
        

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
        real(rkind), dimension(:,:,:), pointer :: k1, k2
        complex(rkind), dimension(:,:,:), pointer :: cbuffz1, cbuffz2

        ! Call this subroutine only after all uhat and u fields are the latest
        ! ones

        k1 => this%spectC%k1
        k2 => this%spectC%k2
        cbuffz1 => this%cbuffzC(:,:,:,1)
        cbuffz2 => this%cbuffzC(:,:,:,2)

        ! Compute omega_x hat
        call transpose_y_to_z(this%vhat,cbuffz1,this%sp_gpC)
        call this%Ops%ddz_C2C(cbuffz1,cbuffz2,topBC_v,botBC_v)
        call transpose_z_to_y(cbuffz2,this%oxhat,this%sp_gpC)
        this%oxhat = imi*k2*this%whatC - this%oxhat 

        ! Compute omega_y hat
        call transpose_y_to_z(this%uhat,cbuffz1,this%sp_gpC)
        call this%Ops%ddz_C2C(cbuffz1,cbuffz2,topBC_u,botBC_u)
        call transpose_z_to_y(cbuffz2,this%oyhat,this%sp_gpC)
        this%oyhat = this%oyhat - imi*k1*this%whatC 

        ! Compute omega_z hat
        this%ozhat = k1*this%vhat 
        this%ozhat = this%ozhat - k2*this%uhat
        this%ozhat = imi*this%ozhat 
   
        ! Compute ox, oy, oz
        call this%spectC%ifft(this%oxhat,this%ox)
        call this%spectC%ifft(this%oyhat,this%oy)
        call this%spectC%ifft(this%ozhat,this%oz)
        
        nullify(k1, k2, cbuffz1, cbuffz2)

    end subroutine


    subroutine addNonLinearTerm_Rot(this)
        class(igrid), intent(inout), target :: this
        real(rkind), dimension(:,:,:), pointer :: rtmpx1
        complex(rkind), dimension(:,:,:), pointer :: ctmpz1, ctmpz2

        
        ! Assume that the rhs vectors haven't been populated

        rtmpx1 => this%rbuffxC(:,:,:,1)
        ctmpz1 => this%cbuffzC(:,:,:,1)
        ctmpz2 => this%cbuffzE(:,:,:,1)
        
        ! x equation
        rtmpx1 = this%oz*this%v
        rtmpx1 = rtmpx1 - this%oy*this%wC
        call this%spectC%fft(rtmpx1,this%u_rhs)

        ! y equation
        rtmpx1 = this%ox*this%wC
        rtmpx1 = rtmpx1 - this%oz*this%u
        call this%spectC%fft(rtmpx1,this%v_rhs)

        ! z equation 
        rtmpx1 = this%oy*this%u
        rtmpx1 = rtmpx1 - this%ox*this%v
        call this%spectC%fft(rtmpx1,this%wC_rhs)
        
        ! Interpolate w_rhs using wC_rhs
        call transpose_y_to_z(this%wC_rhs,ctmpz1, this%sp_gpC)
        call this%Ops%InterpZ_Cell2Edge(ctmpz1,ctmpz2,zeroC,zeroC)   
        call transpose_z_to_y(ctmpz2,this%w_rhs,this%sp_gpE)

        ! Done 
        nullify(rtmpx1) 
        nullify(ctmpz1) 
        nullify(ctmpz2) 
    end subroutine

    subroutine addNonLinearTerm_SkewSymm(this, useCnsrv)
        class(igrid), intent(inout), target :: this
        real(rkind), dimension(:,:,:), pointer :: rtmpx1
        complex(rkind), dimension(:,:,:), pointer :: ctmpz1, ctmpz2!, ctmpz3
        complex(rkind), dimension(:,:,:), pointer :: ctmpy1, ctmpy2
        complex(rkind), dimension(:,:,:), pointer :: dudx, dudy, dudz
        complex(rkind), dimension(:,:,:), pointer :: dvdx, dvdy, dvdz
        complex(rkind), dimension(:,:,:), pointer :: dwdx, dwdy, dwdz
        logical, intent(in) :: useCnsrv
        real(rkind) :: cnst 

        dudx => this%duidxj(:,:,:,1); dudy => this%duidxj(:,:,:,2); dudz => this%duidxj(:,:,:,3); 
        dvdx => this%duidxj(:,:,:,4); dvdy => this%duidxj(:,:,:,5); dvdz => this%duidxj(:,:,:,6); 
        dwdx => this%duidxj(:,:,:,7); dwdy => this%duidxj(:,:,:,8); dwdz => this%duidxj(:,:,:,9); 

        ctmpz1 => this%cbuffzC(:,:,:,1)

        rtmpx1 => this%rbuffxC(:,:,:,1)
        ctmpz1 => this%cbuffzC(:,:,:,1)
        ctmpz2 => this%cbuffzE(:,:,:,1)
        ctmpy1 => this%cbuffyC(:,:,:,1)
        ctmpy2 => this%cbuffyE(:,:,:,1)

        if (useCnsrv) then        
            call this%addNonLinearTerm_Cnsrv()
            cnst = half
            this%u_rhs = half*this%u_rhs
            this%v_rhs = half*this%v_rhs
            this%w_rhs = half*this%w_rhs
        else 
            cnst = one
            this%u_rhs = zero
            this%v_rhs = zero
            this%w_rhs = zero
        end if 


        call this%spectC%ifft(dudx,rtmpx1)
        rtmpx1 = rtmpx1*this%u
        call this%spectC%fft(rtmpx1,ctmpy1)
        this%u_rhs = this%u_rhs - cnst*ctmpy1
        
        call this%spectC%ifft(dudy,rtmpx1)
        rtmpx1 = rtmpx1*this%v
        call this%spectC%fft(rtmpx1,ctmpy1)
        this%u_rhs = this%u_rhs - cnst*ctmpy1

        call this%spectC%ifft(dudz,rtmpx1)
        rtmpx1 = rtmpx1*this%w
        call this%spectC%fft(rtmpx1,ctmpy1)
        this%u_rhs = this%u_rhs - cnst*ctmpy1

        call this%spectC%ifft(dvdx,rtmpx1)
        rtmpx1 = rtmpx1*this%u
        call this%spectC%fft(rtmpx1,ctmpy1)
        this%v_rhs = this%v_rhs - cnst*ctmpy1

        call this%spectC%ifft(dvdy,rtmpx1)
        rtmpx1 = rtmpx1*this%v
        call this%spectC%fft(rtmpx1,ctmpy1)
        this%v_rhs = this%v_rhs - cnst*ctmpy1

        call this%spectC%ifft(dvdz,rtmpx1)
        rtmpx1 = rtmpx1*this%wC
        call this%spectC%fft(rtmpx1,ctmpy1)
        this%v_rhs = this%v_rhs - cnst*ctmpy1

        
        call this%spectC%ifft(dwdx,rtmpx1)
        rtmpx1 = rtmpx1*this%u
        call this%spectC%fft(rtmpx1,ctmpy1)
        call transpose_y_to_z(ctmpy1,ctmpz1,this%sp_gpC)
        call this%Ops%InterpZ_Cell2Edge(ctmpz1,ctmpz2,zeroC,zeroC)
        call transpose_z_to_y(ctmpz2,ctmpy2,this%sp_gpE)  
        this%w_rhs = this%w_rhs - cnst*ctmpy2
    

        call this%spectC%ifft(dwdy,rtmpx1)
        rtmpx1 = rtmpx1*this%v
        call this%spectC%fft(rtmpx1,ctmpy1)
        call transpose_y_to_z(ctmpy1,ctmpz1,this%sp_gpC)
        call this%Ops%InterpZ_Cell2Edge(ctmpz1,ctmpz2,zeroC,zeroC)
        call transpose_z_to_y(ctmpz2,ctmpy2,this%sp_gpE)  
        this%w_rhs = this%w_rhs - cnst*ctmpy2


        call this%spectC%ifft(dwdz,rtmpx1)
        rtmpx1 = rtmpx1*this%wC
        call this%spectC%fft(rtmpx1,ctmpy1)
        call transpose_y_to_z(ctmpy1,ctmpz1,this%sp_gpC)
        call this%Ops%InterpZ_Cell2Edge(ctmpz1,ctmpz2,zeroC,zeroC)
        call transpose_z_to_y(ctmpz2,ctmpy2,this%sp_gpE)  
        this%w_rhs = this%w_rhs - cnst*ctmpy2

        nullify( dudx, dudy, dudz) 
        nullify( dvdx, dvdy, dvdz)
        nullify( dwdx, dwdy, dwdz)
    end subroutine


    subroutine addNonLinearTerm_Cnsrv(this)
        class(igrid), intent(inout), target :: this
        real(rkind), dimension(:,:,:), pointer :: rtmpx1
        complex(rkind), dimension(:,:,:), pointer :: ctmpz1, ctmpz2!, ctmpz3
        complex(rkind), dimension(:,:,:), pointer :: ctmpy1, ctmpy2

        
        ! Assume that the rhs vectors haven't been populated

        rtmpx1 => this%rbuffxC(:,:,:,1)
        ctmpz1 => this%cbuffzC(:,:,:,1)
        ctmpz2 => this%cbuffzE(:,:,:,1)
        !ctmpz3 => this%cbuffzC(:,:,:,2)
        ctmpy1 => this%cbuffyC(:,:,:,1)
        ctmpy2 => this%cbuffyE(:,:,:,1)


        ! uv terms
        rtmpx1 = -this%u*this%v
        call this%spectC%fft(rtmpx1,this%u_rhs)
        this%v_rhs = imi*this%spectC%k1*this%u_rhs
        this%u_rhs = imi*this%spectC%k2*this%u_rhs

        ! uw terms
        rtmpx1 = -this%u*this%wC
        call this%spectC%fft(rtmpx1,ctmpy1)
        call transpose_y_to_z(ctmpy1,ctmpz1,this%sp_gpC)
        call this%Ops%InterpZ_Cell2Edge(ctmpz1,ctmpz2,zeroC,zeroC)
        call transpose_z_to_y(ctmpz2,this%w_rhs,this%sp_gpE)
        this%w_rhs = imi*this%spectE%k1*this%w_rhs
        call this%Ops%ddz_E2C(ctmpz2,ctmpz1)
        call transpose_z_to_y(ctmpz1,ctmpy1,this%sp_gpC)
        this%u_rhs = this%u_rhs + ctmpy1

        ! vw terms
        rtmpx1 = -this%v*this%wC
        call this%spectC%fft(rtmpx1,ctmpy1)
        call transpose_y_to_z(ctmpy1,ctmpz1,this%sp_gpC)
        call this%Ops%InterpZ_Cell2Edge(ctmpz1,ctmpz2,zeroC,zeroC)
        call transpose_z_to_y(ctmpz2,ctmpy2,this%sp_gpE)
        this%w_rhs = this%w_rhs + imi*this%spectE%k2*ctmpy2
        call this%Ops%ddz_E2C(ctmpz2,ctmpz1)
        call transpose_z_to_y(ctmpz1,ctmpy1,this%sp_gpC)
        this%v_rhs = this%v_rhs + ctmpy1
     
        ! uu term 
        rtmpx1 = -this%u*this%u
        call this%spectC%fft(rtmpx1,ctmpy1)
        this%u_rhs = this%u_rhs + this%spectC%k1*ctmpy1
        
        ! vv term
        rtmpx1 = -this%v*this%v
        call this%spectC%fft(rtmpx1,ctmpy1)
        this%v_rhs = this%v_rhs + this%spectC%k2*ctmpy1

        ! ww term 
        rtmpx1 =  -this%wC*this%wC
        call this%spectC%fft(rtmpx1,ctmpy1)
        call transpose_y_to_z(ctmpy1,ctmpz1,this%sp_gpC)
        call this%Ops%ddz_C2E(ctmpz1,ctmpz2,.true., .true.)
        call transpose_z_to_y(ctmpz2,ctmpy2,this%sp_gpE) 
        this%w_rhs = this%w_rhs + ctmpy2

        ! Done 
        nullify(rtmpx1) 
        nullify(ctmpz1) 
        nullify(ctmpz2) 
    end subroutine

    subroutine addCoriolisTerm(this)
        class(igrid), intent(inout) :: this

        ! u equation 
        this%u_rhs = this%u_rhs + this%fCor*this%vhat 

        ! v equation 
        this%v_rhs = this%v_rhs - this%fCor*this%uhat
        if (this%spectC%carryingZeroK) then
            this%v_rhs(this%spectC%ZeroK_i,this%spectC%ZeroK_j,:) =  & 
                    this%v_rhs(this%spectC%ZeroK_i,this%spectC%ZeroK_j,:) & 
                                + this%fCor*this%Gx*this%ny*this%nx
        end if 

        ! w equation 
        ! Do nothing 
    end subroutine  

    subroutine addViscousTerm(this)
        class(igrid), intent(inout), target :: this
        complex(rkind), dimension(:,:,:), pointer :: cytmp1, cztmp1, cztmp2
        complex(rkind), dimension(:,:,:), pointer :: cztmp3, cztmp4, cytmp2

        
        cytmp1 => this%cbuffyC(:,:,:,1)
        cytmp2 => this%cbuffyE(:,:,:,1)
        cztmp1 => this%cbuffzC(:,:,:,1)
        cztmp2 => this%cbuffzC(:,:,:,2)
        cztmp3 => this%cbuffzE(:,:,:,1)
        cztmp4 => this%cbuffzE(:,:,:,2)
        
        ! u equation 
        call transpose_y_to_z(this%uhat, cztmp1,this%sp_gpC)
        call this%Ops%d2dz2_C2C(cztmp1,cztmp2,topBC_u,botBC_u)
        call transpose_z_to_y(cztmp2,cytmp1,this%sp_gpC)
        this%u_rhs = this%u_rhs - this%nu0*this%spectC%kabs_sq *this%uhat 
        this%u_rhs = this%u_rhs + this%nu0*cytmp1
        
        ! v equation 
        call transpose_y_to_z(this%vhat, cztmp1,this%sp_gpC)
        call this%Ops%d2dz2_C2C(cztmp1,cztmp2,topBC_v,botBC_v)
        call transpose_z_to_y(cztmp2,cytmp1,this%sp_gpC)
        this%v_rhs =  this%v_rhs - this%nu0*this%spectC%kabs_sq *this%vhat 
        this%v_rhs =  this%v_rhs + this%nu0*cytmp1 
        
        ! w equation
        call transpose_y_to_z(this%what, cztmp3,this%sp_gpE)
        call this%Ops%d2dz2_E2E(cztmp3,cztmp4,topBC_w,botBC_w)
        call transpose_z_to_y(cztmp4,cytmp2,this%sp_gpC)
        this%w_rhs =  this%w_rhs - this%nu0*this%spectE%kabs_sq *this%what 
        this%w_rhs =  this%w_rhs + this%nu0*cytmp2

        nullify(cytmp1,cytmp2,cztmp1,cztmp2,cztmp3,cztmp4)
    end subroutine 

    subroutine AdamsBashforth(this)
        class(igrid), intent(inout) :: this

        ! Step 1: Non Linear Term 
        select case(AdvectionForm)
        case (1) ! Rotational Form
            call this%AddNonLinearTerm_Rot()
        case (2) ! Conservative Form
            call this%AddNonLinearTerm_Cnsrv()
        case (3) ! Skew Symmetric Form
            call this%AddNonLinearTerm_SkewSymm(.true.)
        case (4) ! Skew Symmetric Form
            call this%AddNonLinearTerm_SkewSymm(.false.)
        end select  

        ! Step 2: Coriolis Term
        if (this%useCoriolis) then
            call this%AddCoriolisTerm()
        end if 

        ! Step 3: Viscous Term
        call this%AddViscousTerm()

        ! Step 4: Time Step 
        if (this%step == 0) then
            this%uhat = this%uhat + this%dt*this%u_rhs
            this%vhat = this%vhat + this%dt*this%v_rhs
            this%what = this%what + this%dt*this%w_rhs
        else
            this%uhat = this%uhat + this%dtby2*(three*this%u_rhs - this%u_Orhs) 
            this%vhat = this%vhat + this%dtby2*(three*this%v_rhs - this%v_Orhs) 
            this%what = this%what + this%dtby2*(three*this%w_rhs - this%w_Orhs) 
        end if 

        ! Step 5: Dealias 
        this%uhat = this%uhat*this%spectC%Gdealias 
        this%vhat = this%vhat*this%spectC%Gdealias 
        this%what = this%what*this%spectE%Gdealias 

        ! Step 6: Pressure projection
        call this%poiss%PressureProjNP(this%uhat,this%vhat,this%what)
        call this%poiss%DivergenceCheck(this%uhat, this%vhat, this%what, this%divergence) 

        ! Step 7: Take it back to physical fields
        call this%spectC%ifft(this%uhat,this%u)
        call this%spectC%ifft(this%vhat,this%v)
        call this%spectE%ifft(this%what,this%w)

        ! STEP 8: Interpolate the cell center values of w
        call this%interp_w_to_wC()
        call this%interp_wHat_to_wHatC()


        ! STEP 9: Compute either vorticity or duidxj 
        select case (AdvectionForm)
        case (1)    
            call this%compute_Vorticity()
        case (3)
            call this%compute_duidxj()
        case (4)
            call this%compute_duidxj()
        end select  

        ! STEP 10: Copy the RHS for using during next time step 
        this%u_Orhs = this%u_rhs
        this%v_Orhs = this%v_rhs
        this%w_Orhs = this%w_rhs

    end subroutine

    subroutine compute_duidxj(this)
        class(igrid), intent(inout), target :: this
        complex(rkind), dimension(:,:,:), pointer :: ctmpz1, ctmpz2
        complex(rkind), dimension(:,:,:), pointer :: ctmpz3, ctmpz4
        complex(rkind), dimension(:,:,:), pointer :: ctmpy1, ctmpy2
        complex(rkind), dimension(:,:,:), pointer :: dudx, dudy, dudz
        complex(rkind), dimension(:,:,:), pointer :: dvdx, dvdy, dvdz
        complex(rkind), dimension(:,:,:), pointer :: dwdx, dwdy, dwdz

        dudx => this%duidxj(:,:,:,1); dudy => this%duidxj(:,:,:,2); dudz => this%duidxj(:,:,:,3); 
        dvdx => this%duidxj(:,:,:,4); dvdy => this%duidxj(:,:,:,5); dvdz => this%duidxj(:,:,:,6); 
        dwdx => this%duidxj(:,:,:,7); dwdy => this%duidxj(:,:,:,8); dwdz => this%duidxj(:,:,:,9); 

        ctmpz1 => this%cbuffzC(:,:,:,1)
        ctmpz2 => this%cbuffzE(:,:,:,1)
        ctmpz3 => this%cbuffzC(:,:,:,2)
        ctmpz4 => this%cbuffzE(:,:,:,2)
        
        ctmpy1 => this%cbuffyC(:,:,:,1)
        ctmpy2 => this%cbuffyE(:,:,:,1)



        dudx = imi*this%spectC%k1*this%uhat
        dudy = imi*this%spectC%k2*this%uhat
        dvdx = imi*this%spectC%k1*this%vhat
        dvdx = imi*this%spectC%k2*this%vhat
        dwdx = imi*this%spectC%k1*this%whatC
        dwdx = imi*this%spectC%k2*this%whatC
       

        call transpose_y_to_z(this%uhat,ctmpz1, this%sp_gpC)
        call this%Ops%ddz_C2C(ctmpz1,ctmpz3,topBC_u,botBC_u) 
        call transpose_z_to_y(ctmpz3,dudz,this%sp_gpC)

        call transpose_y_to_z(this%vhat,ctmpz1, this%sp_gpC)
        call this%Ops%ddz_C2C(ctmpz1,ctmpz3,topBC_v,botBC_v) 
        call transpose_z_to_y(ctmpz3,dvdz,this%sp_gpC)

        call transpose_y_to_z(this%what,ctmpz2, this%sp_gpE)
        call this%Ops%ddz_E2C(ctmpz2,ctmpz3) 
        call transpose_z_to_y(ctmpz3,dwdz,this%sp_gpC)
        nullify( dudx, dudy, dudz) 
        nullify( dvdx, dvdy, dvdz)
        nullify( dwdx, dwdy, dwdz)
    end subroutine
end module 
