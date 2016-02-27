module IncompressibleGridNP
    use kind_parameters, only: rkind, clen
    use constants, only: zero,one,two,three,half 
    use GridMod, only: grid
    use gridtools, only: alloc_buffs, destroy_buffs
    use hooks, only: meshgen, initfields
    use decomp_2d, only: decomp_info, nrank, get_decomp_info, decomp_2d_init, decomp_2d_finalize, &
                    transpose_x_to_y, transpose_y_to_x, transpose_y_to_z, transpose_z_to_y
    use StaggOpsMod, only: staggOps  
    use exits, only: GracefulExit, message
    use spectralMod, only: spectral  
    use PoissonMod, only: poisson
    use mpi 
    use reductions, only: p_maxval

    implicit none

    private
    public :: igrid 

    
    type, extends(grid) :: igrid
        
        type(decomp_info), allocatable :: gpC, gpE
        type(spectral), allocatable :: spectE, spectC
        type(staggOps), allocatable :: Ops

        real(rkind), dimension(:,:,:,:), allocatable :: PfieldsC
        real(rkind), dimension(:,:,:,:), allocatable :: PfieldsE

        complex(rkind), dimension(:,:,:,:), allocatable :: SfieldsC
        complex(rkind), dimension(:,:,:,:), allocatable :: SfieldsE


        type(poisson), allocatable :: poiss
        real(rkind), dimension(:,:,:), allocatable :: divergence

        real(rkind), dimension(:,:,:), pointer :: u, v, w
        complex(rkind), dimension(:,:,:), pointer :: uhat, vhat, what

        real(rkind), dimension(:,:,:,:), allocatable :: rbuffxC, rbuffyC, rbuffzC
        real(rkind), dimension(:,:,:,:), allocatable :: rbuffxE, rbuffyE, rbuffzE
        
        complex(rkind), dimension(:,:,:,:), allocatable :: cbuffxC, cbuffyC, cbuffzC
        complex(rkind), dimension(:,:,:,:), allocatable :: cbuffxE, cbuffyE, cbuffzE

        complex(rkind), dimension(:,:,:,:), allocatable :: duidxj_hat 

        real(rkind) :: nu0, Gx, Gy, Gz

        integer :: runID
        contains
            procedure :: init
            procedure :: destroy
            procedure :: laplacian
            procedure :: gradient 
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
        integer :: t_restartDump = 99999
        logical :: ViscConsrv = .TRUE. 
        logical :: use2DecompFFT = .TRUE.
        logical :: useForcing = .FALSE.
        real(rkind) :: KFmax = 2._rkind  
        namelist /INPUT/       nx, ny, nz, tstop, dt, CFL, nsteps, &
                                              inputdir, outputdir, &
                                  periodicx, periodicy, periodicz, &
                         derivative_x, derivative_y, derivative_z, &
                                     filter_x, filter_y, filter_z, &
                                                       prow, pcol, &
                                             ViscConsrv, SkewSymm, &
                                        t_restartDump, t_dataDump
        namelist /IINPUT/  nu, useSGS, runID , useForcing, KFmax, u_g, use2DecompFFT

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
        call allocate(this%mesh(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3),3))
        call meshgen(this%gpC, this%dx, this%dy, &
            this%dz, this%mesh) ! <-- this procedure is part of user defined HOOKS
        
        ! STEP 4: ALLOCATE/INITIALIZE THE SPECTRAL DERIVED TYPES
        allocate(this%spectC)
        call this%spectC%init("x", nx, ny, nz, this%dx, this%dy, this%dz, &
                this%derivative_x, this%filter_x, 2 , .false.)
        allocate(this%spectE)
        call this%spectE%init("x", nx, ny, nz+1, this%dx, this%dy, this%dz, &
                this%derivative_x, this%filter_x, 2 , .false.)

        ! STEP 5: ALLOCATE/INITIALIZE THE OPERATORS DERIVED TYPE
        call this%Ops%init(this%gpC,this%gpE,0,this%dx,this%dy,this%dz,this%spectC%spectdecomp,this%spectE%spectdecomp)
        

        ! STEP 6: ALLOCATE MEMORY FOR FIELD ARRAYS
        allocate(this%PfieldsC(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3),2))
        allocate(this%divergence(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)))
        allocate(this%PfieldsE(this%gpE%xsz(1),this%gpE%xsz(2),this%gpE%xsz(3),1))
        call this%spectC%alloc_r2c_out(this%SfieldsC,2)
        call this%spectC%alloc_r2c_out(this%duidxj_hat,9)
        call this%spectE%alloc_r2c_out(this%SfieldsE,1)
        
        this%u => this%PfieldsC(:,:,:,1) 
        this%v => this%PfieldsC(:,:,:,2) 
        this%w => this%PfieldsE(:,:,:,1) 
        
        this%uhat => this%SfieldsC(:,:,:,1)
        this%vhat => this%SfieldsC(:,:,:,2)
        this%what => this%SfieldsE(:,:,:,1)

        ! STEP 6: ALLOCATE/INITIALIZE THE POISSON DERIVED TYPE 
        call this%poiss%init(this%spectC,.false.,this%dx,this%dy,this%dz,this%Ops,this%spectE)  
        
        ! STEP 7: INITIALIZE THE FIELDS 
        call initfields_stagg(this%gpC, this%gpE, this%dx, this%dy, this%dz, &
            inputdir, this%mesh, this%PfieldsC, this%PfieldsE)! <-- this procedure is part of user defined HOOKS

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

        ! Check divergence 
        call this%poiss%DivergenceCheck(this%uhat, this%vhat, this%what, this%divergence)
        print*, p_maxval(this%divergence)

        ! STEP 8: ALLOCATE STORAGE FOR BUFFERS AND DUIDXJHAT


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


    end subroutine
end module 
