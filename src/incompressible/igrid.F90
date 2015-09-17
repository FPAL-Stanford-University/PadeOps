module IncompressibleGrid
    use kind_parameters, only: rkind, clen
    use constants, only: zero,one,two,three,half 
    use GridMod, only: grid
    use gridtools, only: alloc_buffs, destroy_buffs
    use hooks, only: meshgen, initfields
    use decomp_2d, only: decomp_info, nrank, get_decomp_info, decomp_2d_init, decomp_2d_finalize, &
                    transpose_x_to_y, transpose_y_to_x, transpose_y_to_z, transpose_z_to_y
    use DerivativesMod,  only: derivatives
  
    use exits, only: GracefulExit, message
    use spectralMod, only: spectral  
    use PoissonMod, only: poisson

    use mpi 

    implicit none

    integer, parameter :: u_index      = 1
    integer, parameter :: v_index      = 2
    integer, parameter :: w_index      = 3
    integer, parameter :: nut_index    = 4

    integer :: numRealBuffs = 2
    integer :: numCmplxBuffs = 2
    character(len=1) :: base_pencil = "x" ! Physical space pencil
    integer :: numFields = 4              ! Number of physical fields to store
    integer :: numT_RHS = 2               ! Number of RHS to store (time steps)

    logical :: fixOddball = .TRUE. 

    type, extends(grid) :: igrid
        ! Spectral realization of the fields 
        complex(rkind), dimension(:,:,:,:), allocatable :: Sfields 
        real(rkind),    dimension(:,:,:,:), allocatable :: duidxj

        ! Molecular Viscosity (useful in DNS)
        real(rkind)                                     :: nu0
        
        ! Pointers
        real(rkind),    dimension(:,:,:), pointer :: x
        real(rkind),    dimension(:,:,:), pointer :: y
        real(rkind),    dimension(:,:,:), pointer :: z
        
        real(rkind),    dimension(:,:,:), pointer :: u
        real(rkind),    dimension(:,:,:), pointer :: v
        real(rkind),    dimension(:,:,:), pointer :: w
        real(rkind),    dimension(:,:,:), pointer :: nut
        
        complex(rkind), dimension(:,:,:), pointer :: uhat
        complex(rkind), dimension(:,:,:), pointer :: vhat
        complex(rkind), dimension(:,:,:), pointer :: what

        real(rkind),    dimension(:,:,:), pointer :: Rtmp1
        real(rkind),    dimension(:,:,:), pointer :: Rtmp2
        complex(rkind), dimension(:,:,:), pointer :: Ctmp1
        complex(rkind), dimension(:,:,:), pointer :: Ctmp2
      
        ! Buffers
        real(rkind),    dimension(:,:,:,:), allocatable :: xbuf_r, ybuf_r, zbuf_r
        complex(rkind), dimension(:,:,:,:), allocatable :: xbuf_c, ybuf_c, zbuf_c
        
        ! Types specific to Incompressible Flow 
        type(spectral), allocatable :: spect
        type(poisson), allocatable :: poiss

        ! LES mode
        logical :: useSGS

        ! Vertical direction BC
        logical :: isZperiodic

        ! RHS array
        complex(rkind), dimension(:,:,:,:,:), allocatable :: rhs_c 
        real(rkind),    dimension(:,:,:,:,:), allocatable :: rhs_r
        
        complex(rkind), dimension(:,:,:,:), pointer :: rhs, rhs_old

        integer :: runID 


        character(len=clen) :: inputdir 
        logical :: use2DecompFFT
        contains
            procedure :: init
            procedure :: destroy
            procedure :: laplacian
            procedure :: gradient
            procedure,private :: CnsrvFrm 
            procedure,private :: NonCnsrvFrm
            procedure,private :: getRHS
            procedure,private :: getDT
            procedure,private :: addVisc
            procedure,private :: get_duidxj
            procedure :: AdamsBashforth
            procedure :: printDivergence

    end type

contains
    subroutine init(this, inputfile )
        class(igrid),target, intent(inout) :: this
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
        integer :: runID = 0
        integer :: t_dataDump = 99999
        integer :: t_restartDump = 99999
        logical :: ViscConsrv = .TRUE. 
        logical :: use2DecompFFT = .TRUE. 
        namelist /INPUT/       nx, ny, nz, tstop, dt, CFL, nsteps, &
                                              inputdir, outputdir, &
                                  periodicx, periodicy, periodicz, &
                         derivative_x, derivative_y, derivative_z, &
                                     filter_x, filter_y, filter_z, &
                                                       prow, pcol, &
                                             ViscConsrv, SkewSymm, &
                                        t_restartDump, t_dataDump
        namelist /IINPUT/  nu, useSGS, runID , use2DecompFFT

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
        this%use2DecompFFT = use2DecompFFT
        this%tsim = zero
        this%tstop = tstop
        this%CFL = CFL

        this%dt = dt 
        this%useSGS = useSGS
        this%step = 0
        this%nsteps = nsteps
        this%t_dataDump = t_dataDump
        this%t_restartDump = t_restartDump

        this%derivative_x = derivative_x
        this%derivative_y = derivative_y
        this%derivative_z = derivative_z

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

        if (.not. periodicx) then
            call GracefulExit("Currently only Periodic BC is supported in x direction",102)
        end if 

        if (.not. periodicy) then
            call GracefulExit("Currently only Periodic BC is supported in y direction",102)
        end if 

        this%isZperiodic = periodicz


        ! STEP 2: INITIALIZE MAIN (PHYSICAL SPACE) DECOMPOSITION
        allocate (this%decomp) 
        call decomp_2d_init(nx, ny, nz, prow, pcol)
        call get_decomp_info(this%decomp)
       

        ! STEP 3: CREATE THE MESH (PHYSICAL SPACE) 
        if ( allocated(this%mesh) ) deallocate(this%mesh) 
        call alloc_buffs(this%mesh,3,base_pencil,this%decomp)
        
        call meshgen(this%decomp, this%dx, this%dy, &
            this%dz, this%mesh) ! <-- this procedure is part of user defined HOOKS
 

        ! STEP 4: INITIALIZE THE "SPECTRAL" DERIVED TYPE 
        if (this%isZperiodic) then 
            allocate(this%spect)
            call this%spect%init(base_pencil, nx, ny, nz, this%dx, this%dy, this%dz, &
                    this%derivative_x, this%filter_x, 3, fixOddball,this%use2DecompFFT, this%ViscConsrv)
        else
            call GracefulExit("CODE INCOMPLETE: Code for non-periodic Z is not &
            & yet complete",102)
        end if 

        ! STEP 5: ALLOCATE FIELDS (BOTH PHYSICAL AND SPECTRAL) 
        if ( allocated(this%fields) ) deallocate(this%fields) 
        call alloc_buffs(this%fields,numfields,base_pencil,this%decomp)
        call alloc_buffs(this%duidxj,9,base_pencil,this%decomp)
        call this%spect%alloc_r2c_out(this%Sfields,3)
        
        this%fields  = zero
        this%Sfields = zero
        this%u  => this%fields(:,:,:,u_index)
        this%v  => this%fields(:,:,:,v_index)
        this%w  => this%fields(:,:,:,w_index)
        this%nut => this%fields(:,:,:,nut_index)
        
        this%uhat => this%Sfields(:,:,:,1)
        this%vhat => this%Sfields(:,:,:,2)
        this%what => this%Sfields(:,:,:,3)


        ! STEP 6: INITIALIZE THE POISSON DERIVED TYPE
        allocate(this%poiss)
        call this%poiss%init(this%spect)


        ! STEP 7: INITIALIZE FIELDS (BOTH PHYSICAL AND SPECTRAL)
        call initfields(this%decomp, this%dx, this%dy, this%dz, &
            inputdir, this%mesh, this%fields)! <-- this procedure is part of user defined HOOKS

        this%nu0 = nu ! <-- initialize nu from the input file 
        
        call this%spect%fft(this%u,this%uhat) ! <-- Convert Phys-to-Spect
        call this%spect%fft(this%v,this%vhat)
        call this%spect%fft(this%w,this%what)

        this%uhat = this%uhat*this%spect%Gdealias ! <-- Dealias/Filter initialized field
        this%vhat = this%vhat*this%spect%Gdealias
        this%what = this%what*this%spect%Gdealias
   
        call this%poiss%PressureProj(this%Sfields,this%spect) !<-- Initial Pressure projection

        call this%spect%ifft(this%uhat,this%u) ! <-- Convert back Spect-to-Phys
        call this%spect%ifft(this%vhat,this%v)
        call this%spect%ifft(this%what,this%w)


        ! STEP 8: INITIALIZE DERIVATIVE DERIVED TYPE (USED FOR POST-PROCESSING) 
        allocate (this%der)
        call this%der%init(                           this%decomp, &
                           this%dx,       this%dy,        this%dz, &
                         periodicx,     periodicy,      periodicz, &
                      derivative_x,  derivative_y,   derivative_z  )      


        ! STEP 9: ALLOCATE BUFFERS AND PHYSICAL SPACE FIELD SIZES 
        call alloc_buffs(this%xbuf_r,numRealBuffs,"x",this%decomp)
        call alloc_buffs(this%ybuf_r,numRealBuffs,"y",this%decomp)
        call alloc_buffs(this%zbuf_r,numRealBuffs,"z",this%decomp) 
        select case (base_pencil)
        case ("x")
            call this%spect%alloc_r2c_out(this%zbuf_c,numCmplxBuffs) 
            this%Rtmp1 => this%xbuf_r(:,:,:,1)
            this%Rtmp2 => this%xbuf_r(:,:,:,2)
            this%Ctmp1 => this%zbuf_c(:,:,:,1)
            this%Ctmp2 => this%zbuf_c(:,:,:,2)
            this%nxp = size(this%xbuf_r,1)
            this%nyp = size(this%xbuf_r,2)
            this%nzp = size(this%xbuf_r,3)
        case ("z")
            call this%spect%alloc_r2c_out(this%xbuf_c,numCmplxBuffs) 
            this%Rtmp1 => this%zbuf_r(:,:,:,1)
            this%Rtmp2 => this%zbuf_r(:,:,:,2)
            this%Ctmp1 => this%xbuf_c(:,:,:,1)
            this%Ctmp2 => this%xbuf_c(:,:,:,2)
            this%nxp = size(this%zbuf_r,1)
            this%nyp = size(this%zbuf_r,2)
            this%nzp = size(this%zbuf_r,3)
        case ("y")
            call GracefulExit("Y- decomposition for the base pencil is not supported",133)
        end select

        ! STEP 10: ALLOCATE THE RHS ARRAY
        allocate(this%rhs_c(size(this%Sfields,1),size(this%Sfields,2),size(this%Sfields,3),size(this%fields,4),numT_RHS))
        this%rhs     => this%rhs_c(:,:,:,:,1)
        this%rhs_old => this%rhs_c(:,:,:,:,2)

        ! STEP 11: INITIALIZE TIMESTEPS  
        this%tsim = zero
        this%step = 0

        ! STEP 12: FINAL CHECK - DID USER SPECIFY EITHER DT OR CFL?
        if ((this%CFL < 0) .and. (this%dt < 0)) then
            call GracefulExit("Neither CFL not dt were specified in the input file", 123)
        end if 
    end subroutine


    subroutine destroy(this)
        class(igrid), intent(inout) :: this
        call this%der%destroy()
        call this%spect%destroy()
        call this%poiss%destroy() 
        deallocate (this%poiss)
        deallocate (this%spect)
        nullify( this%x )    
        nullify( this%y )
        nullify( this%z )
        nullify( this%u )
        nullify( this%v )
        nullify( this%w )
        nullify( this%nut )
        nullify( this%uhat )
        nullify( this%vhat )
        nullify( this%what )
        nullify( this%Rtmp1)
        nullify( this%Rtmp2)
        nullify( this%Ctmp1)
        nullify( this%Ctmp2)
        if (allocated(this%mesh))   deallocate(this%mesh) 
        if (allocated(this%fields)) deallocate(this%fields) 
        deallocate(this%der)
        deallocate(this%decomp)
        if (allocated(this%xbuf_r)) deallocate(this%xbuf_r)
        if (allocated(this%ybuf_r)) deallocate(this%ybuf_r)
        if (allocated(this%zbuf_r)) deallocate(this%zbuf_r)
        if (allocated(this%xbuf_c)) deallocate(this%xbuf_c)
        if (allocated(this%ybuf_c)) deallocate(this%ybuf_c)
        if (allocated(this%zbuf_c)) deallocate(this%zbuf_c)
    end subroutine

    subroutine gradient(this, f, dfdx, dfdy, dfdz)
        class(igrid),target, intent(inout) :: this
        real(rkind), intent(in),  dimension(this%nxp, this%nyp, this%nzp) :: f
        real(rkind), intent(out), dimension(this%nxp, this%nyp, this%nzp) :: dfdx
        real(rkind), intent(out), dimension(this%nxp, this%nyp, this%nzp) :: dfdy
        real(rkind), intent(out), dimension(this%nxp, this%nyp, this%nzp) :: dfdz

        type(derivatives), pointer :: der
        type(decomp_info), pointer :: decomp
        
        der => this%der
        decomp => this%decomp

         dfdx = f
         dfdy = f
         dfdz = f

    end subroutine 

    subroutine laplacian(this, f, lapf)
        use timer
        class(igrid),target, intent(inout) :: this
        real(rkind), intent(in),  dimension(this%nxp, this%nyp, this%nzp) :: f
        real(rkind), intent(out), dimension(this%nxp, this%nyp, this%nzp) :: lapf
        
        type(derivatives), pointer :: der
        type(decomp_info), pointer :: decomp
        

        der => this%der
        decomp => this%decomp

        lapf = f 
    end subroutine

    subroutine AdamsBashforth(this)
        use constants, only: three 
        class(igrid), target, intent(inout) :: this
        type(spectral), pointer :: spec
        real(rkind) :: dtby2

        spec => this%spect
        
        ! STEP 1: Get time step using CFL
        call this%getDT()
        dtby2 = this%dt/two

        ! STEP 2: Get the RHS (in spectral space)
        call this%getRHS()

        ! STEP 3: Check if 1st time step - if yes, do Euler time step, if no, do Adams-Bash
        if (this%step == 0) then
            this%Sfields = this%Sfields + this%dt*this%rhs
        else
            this%Sfields = this%Sfields + dtby2*(three*this%rhs - this%rhs_old) 
        end if 

        ! STEP 4: Update the old RHS
        this%rhs_old = this%rhs 

        ! STEP 5: Dealias 
        this%uhat = this%uhat*spec%Gdealias
        this%vhat = this%vhat*spec%Gdealias
        this%what = this%what*spec%Gdealias

        ! STEP 6: Pressure projection 
        call this%poiss%PressureProj(this%Sfields,this%spect)

        ! STEP 6: Convert from Spectral -> Physical 
        call spec%ifft(this%uhat,this%u)
        call spec%ifft(this%vhat,this%v)
        call spec%ifft(this%what,this%w)

        nullify( spec)
    end subroutine 

    subroutine getDT(this)
        class(igrid), intent(inout) :: this 

        if (this%CFL < 0) then
            return
        end if 
       
        ! Otherwise Compute dt
        this%dt = 0.00001  
    end subroutine 

    subroutine get_duidxj(this)
        use constants, only: imi
        class(igrid), target, intent(inout) :: this
        complex(rkind), dimension(:,:,:), pointer :: Ctmp1, Ctmp2
        type(spectral), pointer :: spec
        real(rkind), dimension(:,:,:), pointer :: dudx, dudy, dudz, dvdx, dvdy,&
                                        dvdz, dwdx, dwdy, dwdz

        spec => this%spect
        ctmp1 => this%zbuf_c(:,:,:,1) ; ctmp2 => this%zbuf_c(:,:,:,2)
        dudx => this%duidxj(:,:,:,1); dudy => this%duidxj(:,:,:,2); dudz => this%duidxj(:,:,:,3)
        dvdx => this%duidxj(:,:,:,4); dvdy => this%duidxj(:,:,:,5); dvdz => this%duidxj(:,:,:,6)
        dwdx => this%duidxj(:,:,:,7); dwdy => this%duidxj(:,:,:,8); dwdz => this%duidxj(:,:,:,9)


        ! STEP 1: Compute the U derivatives
        ctmp1 = imi*this%uhat

        ctmp2 = spec%k1*ctmp1
        call spec%ifft(ctmp2,dudx)

        ctmp2 = spec%k2*ctmp1
        call spec%ifft(ctmp2,dudy)

        ctmp2 = spec%k3*ctmp1
        call spec%ifft(ctmp2,dudz)


        ! STEP 2: Compute the V derivatives
        ctmp1 = imi*this%vhat

        ctmp2 = spec%k1*ctmp1
        call spec%ifft(ctmp2,dvdx)

        ctmp2 = spec%k2*ctmp1
        call spec%ifft(ctmp2,dvdy)

        ctmp2 = spec%k3*ctmp1
        call spec%ifft(ctmp2,dvdz)

        ! STEP 3: Compwte the W derivatives
        ctmp1 = imi*this%what

        ctmp2 = spec%k1*ctmp1
        call spec%ifft(ctmp2,dwdx)

        ctmp2 = spec%k2*ctmp1
        call spec%ifft(ctmp2,dwdy)

        ctmp2 = spec%k3*ctmp1
        call spec%ifft(ctmp2,dwdz)

        nullify( spec, ctmp1, ctmp2)
        nullify( dudx, dudy, dudz, dvdx, dvdy,&
                                        dvdz, dwdx, dwdy, dwdz)
    end subroutine

    subroutine getRHS(this)
        class(igrid), target, intent(inout) :: this

        if ((this%useSGS) .or. (this%SkewSymm)) then 
            call this%get_duidxj() 
        end if 

        call this%CnsrvFrm()
        
        if (this%SkewSymm) then
            call this%NonCnsrvFrm()
        end if 
        
        call this%addVisc()

    end subroutine


    subroutine NonCnsrvFrm(this)
        use constants, only: half 
        class(igrid),target, intent(inout) :: this
        type(spectral), pointer :: spec
        real(rkind), dimension(:,:,:), pointer :: tmp1
        complex(rkind), dimension(:,:,:), pointer :: ctmp1
        real(rkind), dimension(:,:,:), pointer :: dudx, dudy, dudz, dvdx, dvdy,&
                                        dvdz, dwdx, dwdy, dwdz


        spec => this%spect
        tmp1 => this%xbuf_r(:,:,:,1)
        ctmp1 => this%zbuf_c(:,:,:,1)
        dudx => this%duidxj(:,:,:,1); dudy => this%duidxj(:,:,:,2); dudz => this%duidxj(:,:,:,3)
        dvdx => this%duidxj(:,:,:,4); dvdy => this%duidxj(:,:,:,5); dvdz => this%duidxj(:,:,:,6)
        dwdx => this%duidxj(:,:,:,7); dwdy => this%duidxj(:,:,:,8); dwdz => this%duidxj(:,:,:,9)


        ! STEP 1: Halve the existing RHS
        this%rhs = half*this%rhs

        ! Step 2: Add the contribution to u equation
        tmp1 = this%u*dudx
        tmp1 = tmp1 + this%v*dudy
        tmp1 = tmp1 + this%w*dudz
        call spec%fft(tmp1,ctmp1)
        this%rhs(:,:,:,1) = this%rhs(:,:,:,1) + half*ctmp1


        ! Step 3: Add the contribution to v equation
        tmp1 = this%u*dvdx
        tmp1 = tmp1 + this%v*dvdy
        tmp1 = tmp1 + this%w*dvdz
        call spec%fft(tmp1,ctmp1)
        this%rhs(:,:,:,2) = this%rhs(:,:,:,2) + half*ctmp1

        ! Step 3: Add the contribution to w equation
        tmp1 = this%u*dwdx
        tmp1 = tmp1 + this%v*dwdy
        tmp1 = tmp1 + this%w*dwdz
        call spec%fft(tmp1,ctmp1)
        this%rhs(:,:,:,3) = this%rhs(:,:,:,3) + half*ctmp1

        nullify(spec, ctmp1, tmp1)
        nullify( dudx, dudy, dudz, dvdx, dvdy,&
                                        dvdz, dwdx, dwdy, dwdz)
    end subroutine

    subroutine CnsrvFrm(this)
        use constants, only: imi 
        class(igrid), target, intent(inout) :: this
        type(spectral), pointer :: spec
        real(rkind), dimension(:,:,:), pointer :: tmp1
        complex(rkind), dimension(:,:,:), pointer :: ctmp1

        spec => this%spect
        tmp1 => this%xbuf_r(:,:,:,1)
        ctmp1 => this%zbuf_c(:,:,:,1)


        ! STEP 1: Compute uu derivatives and add
        tmp1 = this%u*this%u
        call spec%fft(tmp1,ctmp1) 
        this%rhs(:,:,:,1) = -spec%k1*ctmp1 

        ! STEP 2: Compute uv derivatives and add
        tmp1 = this%u*this%v
        call spec%fft(tmp1,ctmp1) 
        this%rhs(:,:,:,1) = this%rhs(:,:,:,1) - spec%k2*ctmp1 
        this%rhs(:,:,:,2) = -spec%k1*ctmp1 

        ! STEP 3: Compute uw derivatives and add
        tmp1 = this%u*this%w
        call spec%fft(tmp1,ctmp1) 
        this%rhs(:,:,:,1) = this%rhs(:,:,:,1) - spec%k3*ctmp1 
        this%rhs(:,:,:,3) = -spec%k1*ctmp1 

        ! STEP 4: Compute vv derivatives and add
        tmp1 = this%v*this%v
        call spec%fft(tmp1,ctmp1) 
        this%rhs(:,:,:,2) = this%rhs(:,:,:,2) - spec%k2*ctmp1 

        ! STEP 5: Compute vw derivatives and add
        tmp1 = this%v*this%w
        call spec%fft(tmp1,ctmp1) 
        this%rhs(:,:,:,2) = this%rhs(:,:,:,2) - spec%k3*ctmp1 
        this%rhs(:,:,:,3) = this%rhs(:,:,:,3) - spec%k2*ctmp1 

        ! STEP 6: Compute ww derivatives and add
        tmp1 = this%w*this%w
        call spec%fft(tmp1,ctmp1) 
        this%rhs(:,:,:,3) = this%rhs(:,:,:,3) - spec%k3*ctmp1 
   
        ! STEP 7: Multiply RHS by 1i
        this%rhs = imi*this%rhs 

        nullify(spec, ctmp1, tmp1)
    end subroutine

    subroutine addVisc(this)
        class(igrid),target, intent(inout) :: this
        type(spectral), pointer :: spec
        complex(rkind), dimension(:,:,:), pointer :: ctmp1, ctmp2

        ! if LES the update nut and that part of RHS here
        
        ctmp1 => this%zbuf_c(:,:,:,1)
        ctmp2 => this%zbuf_c(:,:,:,2)
        spec => this%spect

        ! STEP 1: Compute the laplacian multiplier
        ctmp1 =  this%nu0 * ( spec%k1_der2 + spec%k2_der2  + spec%k3_der2 )

        ! STEP 2: Update the u equation
        ctmp2 = ctmp1*this%uhat 
        this%rhs(:,:,:,1) = this%rhs(:,:,:,1) - ctmp2 
        
        ! STEP 3: Update the v equation
        ctmp2 = ctmp1*this%vhat 
        this%rhs(:,:,:,2) = this%rhs(:,:,:,2) - ctmp2 
        
        ! STEP 4: Update the u equation
        ctmp2 = ctmp1*this%what 
        this%rhs(:,:,:,3) = this%rhs(:,:,:,3) - ctmp2 

        nullify(spec, ctmp1, ctmp2)
    end subroutine 

    subroutine printDivergence(this)
        use reductions, only: p_maxval
        use constants, only: imi
        class(igrid),target, intent(inout) :: this
        type(derivatives), pointer :: der
        
        real(rkind), dimension(:,:,:), pointer :: xtmp1, xtmp2
        real(rkind), dimension(:,:,:), pointer :: ytmp1, ytmp2
        real(rkind), dimension(:,:,:), pointer :: ztmp1, ztmp2
        complex(rkind), dimension(:,:,:), pointer :: Cztmp1
        
        der => this%der 
        xtmp1 => this%xbuf_r(:,:,:,1) 
        xtmp2 => this%xbuf_r(:,:,:,2) 
        ytmp1 => this%ybuf_r(:,:,:,1) 
        ytmp2 => this%ybuf_r(:,:,:,2) 
        ztmp1 => this%zbuf_r(:,:,:,1) 
        ztmp2 => this%zbuf_r(:,:,:,2) 
       
        ! STEP 1: Get x derivative 
        call der%ddx(this%u,xtmp1) 
        xtmp2 = xtmp1

        ! STEP 2: Get the y derivative and add 
        call transpose_x_to_y(this%v,ytmp1,this%decomp)
        call der%ddy(ytmp1,ytmp2)
        call transpose_y_to_x(ytmp2,xtmp1,this%decomp)
        xtmp2 = xtmp2 + xtmp1

        ! STEP 3: Get the z derivative and add
        call transpose_x_to_y(this%w,ytmp1,this%decomp)
        call transpose_y_to_z(ytmp1,ztmp1,this%decomp)
        call der%ddz(ztmp1,ztmp2)
        call transpose_z_to_y(ztmp2,ytmp2,this%decomp)
        call transpose_y_to_x(ytmp2,xtmp1,this%decomp)
        xtmp2 = xtmp2 + xtmp1
       
        ! STEP 4: Print the maximum value of the divergence  
        call message(2,"Max divergence",p_maxval(xtmp2))


        ! Alternately:
        !Cztmp1 => this%zbuf_c(:,:,:,1)
        !Cztmp1 = imi*(this%uhat*this%spect%k1 + this%vhat*this%spect%k2 + this%what*this%spect%k3)
        !call this%spect%ifft(Cztmp1,xtmp1)
        !call message("Maximum divergence",p_maxval(xtmp1))

        nullify (Cztmp1, der, xtmp1, xtmp2, ytmp1, ytmp2,ztmp1, ztmp2)
    end subroutine

end module 
