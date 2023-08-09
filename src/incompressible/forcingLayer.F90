module forcingLayerMod
  ! This forcing scheme is adapted from that presented in 
  ! Bodart, Julien, J-B. Cazalbou, and Laurent Joly. "Direct numerical
  ! simulation of unsheared turbulence diffusing towards a free-slip or no-slip
  ! surface." Journal of Turbulence 11 (2010): N48.
    use kind_parameters, only: rkind, clen, mpirkind
    use decomp_2d, only: decomp_info, transpose_x_to_y, transpose_y_to_z, &
        transpose_z_to_y, transpose_y_to_x, nrank, nproc
    use decomp_2d_io, only: decomp_2d_read_one
    use spectralMod, only: spectral
    use reductions, only: p_minval, p_maxval, p_sum
    use PadeDerOps, only: Pade6Stagg
    use seedGen, only: initializeLogMap
    use PadePoissonMod, only: Padepoisson 
    use exits, only: message, message_min_max
    use basic_io, only: read_1d_ascii
    use mpi
    use fortran_assert, only: assert
    use constants, only: third, half, rmaxInt, im0
    use interpolatorMod, only: interpolator
    
    implicit none
    private

    real(rkind), dimension(4) :: xsample = [0.d0, 1.d0/3.d0, 2.d0/3.d0, 1.d0]
    real(rkind), dimension(4) :: vsample = [0.d0, 2.d0/3.d0, 2.d0/3.d0, 0.d0]
    integer, parameter :: zbc_bot = 0, zbc_top = 0
    public :: forcingLayer, onThisRank, getStEndIndices

    type :: forcingLayer
        integer :: nblocks ! The number of forcing blocks in each x,y direction
        real(rkind) :: lf  ! Physical size of forcing layer
        real(rkind), dimension(:,:,:), allocatable :: fx, fy, fz
        complex(rkind), dimension(:,:,:), allocatable :: fxhat, fyhat, fzhat
        complex(rkind), dimension(:,:,:), allocatable :: fxhat_old, fyhat_old, fzhat_old
        real(rkind), dimension(:,:,:), allocatable :: phixC, phiyC, phizC
        real(rkind), dimension(:,:,:), allocatable :: phixE, phiyE, phizE
        real(rkind), dimension(:,:,:), allocatable :: dphixC, dphiyC, dphizC
        real(rkind), dimension(:,:,:), allocatable :: dphixE, dphiyE, dphizE
        type(decomp_info), pointer :: gpC, gpE, sp_gpC, sp_gpE
        type(spectral), pointer :: spectC, spectE
        real(rkind), dimension(4) :: xs, ys, zs, vs ! Sample points for spline interpolation
        real(rkind), dimension(:,:,:), pointer :: x, y, z, xE, yE, zE
        real(rkind) :: zgap, zst
        real(rkind), dimension(3) :: periodicXshift, periodicYshift
        type(Pade6Stagg), pointer :: Pade6opZ
        logical :: dumpForce, dumpSplines, checkDivergence
        logical :: randomizeBlockPosition_xy, randomizeBlockPosition_z
        real(rkind) :: seedFact
        integer :: seed, fixedForceType
        real(rkind) :: xshift = 0.d0, yshift = 0.d0
        real(rkind), dimension(:,:), allocatable :: zshift
        type(padePoisson), pointer :: poiss
        integer, dimension(:,:), allocatable :: forceType
        real(rkind) :: Lx, Ly
        real(rkind), dimension(:,:,:), allocatable :: k1C, k1E
        real(rkind), dimension(:,:,:), allocatable :: sp_buffyC, sp_buffyE
        integer :: version
        character(len=clen) :: inputdir
        logical :: initFull = .false.
        real(rkind) :: maxDiv

        ! Required for forcing control strategies:
        integer :: controlMethod
        real(rkind) :: ampFact, tgtDissipation, tgtKE, integralTime, gain
        integer :: kst, ken

        contains
          procedure, private :: init_full
          procedure, private :: init_lightweight
          generic            :: init => init_full, init_lightweight
          procedure :: reinit
          procedure :: destroy
          procedure :: updateRHS
          procedure :: updateForce
          procedure :: interpC2E
          procedure :: interpE2C
          procedure :: updateSeed
          procedure :: divergenceCheck
          procedure :: readRestartFiles
    end type

    contains

      subroutine init_lightweight(this,gpC,gpE,TID,runID,restartSim,inputdir)
          class(forcingLayer), intent(inout) :: this
          type(decomp_info), intent(inout), target :: gpC, gpE
          integer, intent(in) :: TID, runID
          logical, intent(in) :: restartSim
          character(len=*), intent(in) :: inputdir

          this%gpC => gpC
          this%gpE => gpE
          this%inputdir = inputdir
          
          allocate(this%fx( this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3)))
          allocate(this%fy( this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3)))
          allocate(this%fz( this%gpE%xsz(1), this%gpE%xsz(2), this%gpE%xsz(3)))
          
          if (restartSim) then
              call this%readRestartFiles(runID,TID,this%fx,this%fy,this%fz,&
                this%gpC,this%gpE)
          end if

          this%initFull = .false.
      end subroutine

      subroutine init_full(this,inputfile,Lx,Ly,mesh,xE,yE,zE,gpC,gpE,spectC,spectE,Pade6opZ,PadePoiss,&
          runID,TID,inputdir,rbuffxC,rbuffxE,cbuffxC,cbuffxE,cbuffyC,cbuffyE,cbuffzC,cbuffzE,restartSim)
          class(forcingLayer), intent(inout) :: this
          character(len=*), intent(in) :: inputfile, inputdir
          real(rkind), intent(in) :: Lx, Ly
          real(rkind), dimension(:,:,:,:), intent(in), target :: mesh ! Cell-centered mesh (x,y,z)
          real(rkind), dimension(:,:,:), intent(in), target :: xE, yE, zE ! Edge mesh (z only)
          class(decomp_info), intent(in), target :: gpC, gpE
          class(spectral), intent(in), target :: spectC, spectE
          class(Pade6Stagg), intent(in), target :: Pade6opZ
          class(padePoisson), intent(in), target :: PadePoiss
          logical, intent(in) :: restartSim
          integer, intent(in) :: runID, TID
          real(rkind), dimension(:,:,:,:), intent(inout) :: rbuffxC, rbuffxE
          complex(rkind), dimension(:,:,:,:), intent(inout) :: cbuffyC, cbuffyE, cbuffzC, cbuffzE
          complex(rkind), dimension(:,:,:), intent(inout) :: cbuffxC, cbuffxE
          integer :: nblocks = 1
          integer :: ioUnit, ierr, fixedForceType = 0
          real(rkind) :: zgap = 0
          logical :: dumpForce = .false., dumpSplines = .false., checkDivergence = .false.
          logical :: randomizeBlockPosition_xy = .true.
          logical :: randomizeBlockPosition_z = .true.
          real(rkind), dimension(:,:,:), allocatable :: tmp
          integer :: version = 1, controlMethod = 1
          real(rkind) :: tgtDissipation = -1.d0, tgtKE = -1.d0, integralTime = -1.d0, gain = 50.d0
          real(rkind), dimension(:), allocatable :: z1D
          real(rkind) :: Lf

          namelist /localizedForceLayer/ nblocks, zgap, dumpForce, &
            dumpSplines, checkDivergence, fixedForceType, randomizeBlockPosition_xy, &
            randomizeBlockPosition_z, version, tgtDissipation, tgtKE, gain, &
            controlMethod
          
          ioUnit = 123
          open(unit=ioUnit, file=trim(inputfile), form='FORMATTED', iostat=ierr)
          read(unit=ioUnit, NML=localizedForceLayer)
          close(ioUnit)

          ! Input file parameters
          !     nblocks --> Number of forcing blocks in each x, y direction
          !               we can randomly shift the forcing blocks in x and y without worrying
          !               about periodicity
          !     zgap --> Allows random shifts in z

          call message(0,'=========== Initializing forcing layer ============')
          this%nblocks  = nblocks 
          this%zgap     = zgap
          this%lf       = Lx/real(nblocks,rkind)
          this%Lx       = Lx
          this%Ly       = Ly

          this%gpC => gpC 
          this%gpE => gpE
          this%spectC => spectC
          this%spectE => spectE
          this%sp_gpC => spectC%spectdecomp
          this%sp_gpE => spectE%spectdecomp
          this%Pade6opZ => Pade6opZ
          this%poiss => PadePoiss

          this%xs = xsample*this%lf
          this%ys = this%xs
          this%zs = this%xs
          this%vs = vsample
          
          this%x => mesh(:,:,:,1)
          this%y => mesh(:,:,:,2)
          this%z => mesh(:,:,:,3)
          this%xE => xE
          this%yE => yE
          this%zE => zE

          this%dumpForce = dumpForce
          this%dumpSplines = dumpSplines
          this%checkDivergence = checkDivergence
          this%randomizeBlockPosition_xy = randomizeBlockPosition_xy
          this%randomizeBlockPosition_z = randomizeBlockPosition_z

          this%fixedForceType = fixedForceType
          this%version = version ! 1: Same fundamental force types as Bodart et al.
                                 ! 2: Force types have zero mean locally

          this%inputdir = inputdir

          ! The following is used to scale the forcing term such that a target
          ! stationary state is achieved
!print*, "conrolMethod =", controlMethod
!call MPI_Barrier(MPI_COMM_WORLD,ierr)
          call assert(controlMethod < 4,'Must select appropriate control method -- forcingLayer.F90')
          call assert(controlMethod > 0,'Must select appropriate control method -- forcingLayer.F90')
          call assert(tgtDissipation > 0.d0,'Must set tgtDissipation to a postivive value -- forcingLayer.F90')
          if (controlMethod == 2)  then
              call assert(tgtKE > 0.d0,'Must set tgtKE to a postivive value -- forcingLayer.F90')
          end if
          this%controlMethod  = controlMethod
          this%tgtDissipation = tgtDissipation
          this%tgtKE          = tgtKE
          this%integralTime   = tgtKE/tgtDissipation
          this%gain           = gain
        
          allocate(this%forceType(this%nblocks,this%nblocks))

          allocate(this%fx( this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3)))
          allocate(this%fy( this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3)))
          allocate(this%fz( this%gpE%xsz(1), this%gpE%xsz(2), this%gpE%xsz(3)))

          allocate(this%fxhat(this%spectC%spectDecomp%ysz(1), this%spectC%spectDecomp%ysz(2), this%spectC%spectDecomp%ysz(3)))
          allocate(this%fyhat(this%spectC%spectDecomp%ysz(1), this%spectC%spectDecomp%ysz(2), this%spectC%spectDecomp%ysz(3)))
          allocate(this%fzhat(this%spectE%spectDecomp%ysz(1), this%spectE%spectDecomp%ysz(2), this%spectE%spectDecomp%ysz(3)))
          
          allocate(this%fxhat_old(this%spectC%spectDecomp%ysz(1), this%spectC%spectDecomp%ysz(2), this%spectC%spectDecomp%ysz(3)))
          allocate(this%fyhat_old(this%spectC%spectDecomp%ysz(1), this%spectC%spectDecomp%ysz(2), this%spectC%spectDecomp%ysz(3)))
          allocate(this%fzhat_old(this%spectE%spectDecomp%ysz(1), this%spectE%spectDecomp%ysz(2), this%spectE%spectDecomp%ysz(3)))
          
          allocate(this%phixC(this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3)))
          allocate(this%phiyC(this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3)))
          allocate(this%phizC(this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3)))
          
          allocate(this%phixE(this%gpE%xsz(1), this%gpE%xsz(2), this%gpE%xsz(3)))
          allocate(this%phiyE(this%gpE%xsz(1), this%gpE%xsz(2), this%gpE%xsz(3)))
          allocate(this%phizE(this%gpE%xsz(1), this%gpE%xsz(2), this%gpE%xsz(3)))
          
          allocate(this%dphixC(this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3)))
          allocate(this%dphiyC(this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3)))
          allocate(this%dphizC(this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3)))
          
          allocate(this%dphixE(this%gpE%xsz(1), this%gpE%xsz(2), this%gpE%xsz(3)))
          allocate(this%dphiyE(this%gpE%xsz(1), this%gpE%xsz(2), this%gpE%xsz(3)))
          allocate(this%dphizE(this%gpE%xsz(1), this%gpE%xsz(2), this%gpE%xsz(3)))

          allocate(this%zshift(nblocks,nblocks))
          this%zshift = 0.d0

          allocate(this%k1C(this%sp_gpC%ysz(1),this%sp_gpC%ysz(2),this%sp_gpC%ysz(3)))
          allocate(this%k1E(this%sp_gpE%ysz(1),this%sp_gpE%ysz(2),this%sp_gpE%ysz(3)))
          allocate(this%sp_buffyC(this%sp_gpC%ysz(1),this%sp_gpC%ysz(2),this%sp_gpC%ysz(3)))
          allocate(this%sp_buffyE(this%sp_gpE%ysz(1),this%sp_gpE%ysz(2),this%sp_gpE%ysz(3)))

          allocate(z1D(this%gpC%xsz(3)))
          z1D = this%z(1,1,:)

          this%k1C = spectC%k1
          this%k1E = spectE%k1

          allocate(tmp(this%sp_gpC%xsz(1),this%sp_gpC%xsz(2),this%sp_gpC%xsz(3)))
          call transpose_y_to_x(this%k1C,tmp,this%sp_gpC)
          tmp(spectC%nx_g/2+1,:,:) = 0
          call transpose_x_to_y(tmp,this%k1C,this%sp_gpC)
          deallocate(tmp)
          allocate(tmp(this%sp_gpE%xsz(1),this%sp_gpE%xsz(2),this%sp_gpE%xsz(3)))
          call transpose_y_to_x(this%k1E,tmp,this%sp_gpE)
          tmp(spectE%nx_g/2+1,:,:) = 0
          call transpose_x_to_y(tmp,this%k1E,this%sp_gpE)
          deallocate(tmp)
          
          ! Location of the bottom of the forcing layer
          this%zst = -0.5d0*this%lf

          this%periodicXshift = [0.d0, this%Lx, -this%Lx]
          this%periodicYshift = [0.d0, this%Ly, -this%Ly]

          ! Initialize randomized seeds
          if (restartSim) then
              call this%readRestartFiles(runID,TID,this%fx,this%fy,this%fz,&
                this%gpC,this%gpE)
             
              ! Transform to spectral space
              call this%spectC%fft(this%fx,this%fxhat_old)
              call this%spectC%fft(this%fy,this%fyhat_old)
              call this%spectE%fft(this%fz,this%fzhat_old)

              call fixOddBall(this%fxhat_old,this%spectC,cbuffxC)
              call fixOddBall(this%fyhat_old,this%spectC,cbuffxC)
              call fixOddBall(this%fzhat_old,this%spectE,cbuffxE)
             
              ! Check that the force is divergence free
              this%fxhat = this%fxhat_old 
              this%fyhat = this%fyhat_old 
              this%fzhat = this%fzhat_old 
              
          else
              this%seedFact = initializeLogMap(1)
              this%ampFact = 1.d0
              
              call this%updateForce(rbuffxC, rbuffxE, cbuffxC, cbuffxE, cbuffyC, cbuffyE, cbuffzC, cbuffzE)
              this%fxhat_old = this%fxhat
              this%fyhat_old = this%fyhat
              this%fzhat_old = this%fzhat
              
              call this%spectC%ifft(this%fxhat,this%fx)
              call this%spectC%ifft(this%fyhat,this%fy)
              call this%spectE%ifft(this%fzhat,this%fz)

          end if

          ! Set bounds for computing integrated work, dissipation, etc. for the
          ! forcing control strategies
          Lf = this%lf + this%zgap
          if (onThisRank(z1D,-0.5d0*Lf,0.5d0*Lf)) then
              call getStEndIndices(z1D,-0.5d0*Lf,0.5d0*Lf,this%kst,this%ken)
          else
              this%kst = 2
              this%ken = 1
          end if

          deallocate(z1D)
          
          call message(1,'Forcing layer initialized')
          call message_min_max(1,"Bounds for fx:", this%ampFact*p_minval(minval(this%fx)), &
            this%ampFact*p_maxval(maxval(this%fx)))
          call message_min_max(1,"Bounds for fy:", this%ampFact*p_minval(minval(this%fy)), &
            this%ampFact*p_maxval(maxval(this%fy)))
          call message_min_max(1,"Bounds for fz:", this%ampFact*p_minval(minval(this%fz)), &
            this%ampFact*p_maxval(maxval(this%fz)))
          call this%divergenceCheck(rbuffxC(:,:,:,1),fixDivergence = .false.)

          if (this%maxDiv > 1.d-12) then
              ! The only reason this should happen is if we are restarting from
              ! a courser simulation. In that case fx, fy, and fz are
              ! interpolated from the course grid useing linear interpolation
              ! and so the force will not be divergence free.
              call message(1,"Projecting out divergence")

              call this%divergenceCheck(rbuffxC(:,:,:,1),fixDivergence = .true.)
             
              ! Get physical space force 
              call this%spectC%ifft(this%fxhat,this%fx)
              call this%spectC%ifft(this%fyhat,this%fy)
              call this%spectE%ifft(this%fzhat,this%fz)

              ! Store the old fhats for next time step
              this%fxhat_old = this%fxhat
              this%fyhat_old = this%fyhat
              this%fzhat_old = this%fzhat
            
              ! Print force bounds
              call message_min_max(1,"Bounds for fx:", this%ampFact*p_minval(minval(this%fx)), &
                this%ampFact*p_maxval(maxval(this%fx)))
              call message_min_max(1,"Bounds for fy:", this%ampFact*p_minval(minval(this%fy)), &
                this%ampFact*p_maxval(maxval(this%fy)))
              call message_min_max(1,"Bounds for fz:", this%ampFact*p_minval(minval(this%fz)), &
                this%ampFact*p_maxval(maxval(this%fz)))
          end if
          
          call message(1,'Forcing layer bounds (about forcing mid-plane)',[-0.5d0*(this%zgap + this%lf), &
            0.5d0*(this%zgap + this%lf)])
          
          this%initFull = .true.
      end subroutine
      
      subroutine reinit(this,runID,TID,cbuffxC,cbuffxE)
          class(forcingLayer), intent(inout) :: this
          integer, intent(in) :: runID, TID
          complex(rkind), dimension(:,:,:), intent(inout) :: cbuffxC, cbuffxE
          integer :: ierr

          ! Read force info
          call this%readRestartFiles(runID,TID,this%fx,this%fy,this%fz,&
            this%gpC,this%gpE)

          call this%spectC%fft(this%fx,this%fxhat_old)
          call this%spectC%fft(this%fy,this%fyhat_old)
          call this%spectE%fft(this%fz,this%fzhat_old)

          call fixOddBall(this%fxhat_old,this%spectC,cbuffxC)
          call fixOddBall(this%fyhat_old,this%spectC,cbuffxC)
          call fixOddBall(this%fzhat_old,this%spectE,cbuffxC)

          call message(0,'Forcing layer re-initialized')
          call message(1,'Forcing layer bounds (about forcing mid-plane)',[-0.5d0*(this%zgap + this%lf), &
            0.5d0*(this%zgap + this%lf)])
      end subroutine
          
      subroutine readRestartFiles(this,runID,TID,fx,fy,fz,gpC,gpE)
          class(forcingLayer), intent(inout) :: this
          integer, intent(in) :: runID, TID
          real(rkind), dimension(:,:,:), intent(out) :: fx, fy, fz
          type(decomp_info), intent(in) :: gpC, gpE
          character(len=clen) :: tempname, fname
          real(rkind), dimension(:), allocatable :: tmp
          integer :: ierr

          if (nrank == 0) then
              write(tempname,"(A7,A4,I2.2,A7,I6.6)") "RESTART", "_Run",runID, "_finfo.",TID
              fname = this%inputDir(:len_trim(this%inputDir))//"/"//trim(tempname)
              call read_1d_ascii(tmp,fname)
              this%seedFact = tmp(1)
              this%ampFact = tmp(2)
              deallocate(tmp)
          end if
          call MPI_Barrier(MPI_COMM_WORLD,ierr)
          call MPI_Bcast(this%seedFact,1,mpirkind,0,MPI_COMM_WORLD,ierr)
          call MPI_Bcast(this%ampFact,1,mpirkind,0,MPI_COMM_WORLD,ierr)
          call MPI_Barrier(MPI_COMM_WORLD,ierr)

          write(tempname,"(A7,A4,I2.2,A4,I6.6)") "RESTART", "_Run",runID, "_fx.",TID
          fname = this%inputDir(:len_trim(this%inputDir))//"/"//trim(tempname)
          call decomp_2d_read_one(1,fx, fname, gpC)
          
          write(tempname,"(A7,A4,I2.2,A4,I6.6)") "RESTART", "_Run",runID, "_fy.",TID
          fname = this%inputDir(:len_trim(this%inputDir))//"/"//trim(tempname)
          call decomp_2d_read_one(1,fy, fname, gpC)
          
          write(tempname,"(A7,A4,I2.2,A4,I6.6)") "RESTART", "_Run",runID, "_fz.",TID
          fname = this%inputDir(:len_trim(this%inputDir))//"/"//trim(tempname)
          call decomp_2d_read_one(1,fz, fname, gpE)

      end subroutine

      subroutine destroy(this)
          class(forcingLayer), intent(inout) :: this
         
          if (this%initFull) then
              deallocate(this%fxhat, this%fyhat, this%fzhat)
              deallocate(this%fxhat_old, this%fyhat_old, this%fzhat_old)
              deallocate(this%phixC, this%phiyC, this%phizC)
              deallocate(this%phixE, this%phiyE, this%phizE)
              deallocate(this%dphixC, this%dphiyC, this%dphizC)
              deallocate(this%dphixE, this%dphiyE, this%dphizE)
              deallocate(this%zshift, this%forceType, this%k1C, this%k1E)
              deallocate(this%sp_buffyC, this%sp_buffyE)
          
              nullify(this%sp_gpC, this%sp_gpE)
              nullify(this%spectC, this%spectE)
              nullify(this%x, this%y, this%z, this%xE, this%yE, this%zE)
              nullify(this%Pade6opZ, this%poiss)
          end if 
          
          deallocate(this%fx, this%fy, this%fz)
          nullify(this%gpC, this%gpE)
      end subroutine

      subroutine updateRHS(this,urhs,vrhs,wrhs,u,v,wC,duidxjC,Re,nSGS,&
          cbuffxC,cbuffxE,cbuffyC,cbuffyE,cbuffzC,cbuffzE,&
          rbuffxC,rbuffxE,rbuffyC,rbuffyE,rbuffzC,rbuffzE,newTimeStep,dt)
          ! The force added to the RHS is computed from f = ddt(Bodart force).
          ! Below Bodart force is first computed and stored in fhat, then fhat
          ! is updated via Forward Euler.
          class(forcingLayer), intent(inout), target :: this
          complex(rkind), dimension(this%sp_gpC%ysz(1), this%sp_gpC%ysz(2), this%sp_gpC%ysz(3)), intent(inout) :: urhs, vrhs
          complex(rkind), dimension(this%sp_gpE%ysz(1), this%sp_gpE%ysz(2), this%sp_gpE%ysz(3)), intent(inout) :: wrhs
          real(rkind), dimension(:,:,:), intent(in), target :: u, v, wC
          real(rkind), dimension(:,:,:), intent(in) :: nSGS
          real(rkind), dimension(:,:,:,:), intent(in) :: duidxjC
          real(rkind), intent(in) :: Re
          complex(rkind), dimension(:,:,:), intent(inout) :: cbuffxC, cbuffxE
          complex(rkind), dimension(:,:,:,:), intent(inout) :: cbuffyC, cbuffyE
          real(rkind), dimension(:,:,:,:), intent(inout) :: rbuffxC, rbuffxE, &
            rbuffyC, rbuffyE, rbuffzC, rbuffzE
          logical, intent(in) :: newTimeStep
          real(rkind), intent(in) :: dt
          complex(rkind), dimension(:,:,:,:), intent(inout) :: cbuffzC, cbuffzE
          real(rkind) :: forceWork, KE, eps
          real(rkind), dimension(:,:,:), pointer :: uF, vF, wF, fxF, fyF, fzF, &
            dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz, nSGS_F

          call assert(this%initFull,'Cannot call updateRHS without doing'//&
            ' full initialization -- forcingLayer.F0-')
          if (newTimeStep) then
              ! Link pointers
              associate( uF => u(:,:,this%kst:this%ken), vF => v(:,:,this%kst:this%ken), &
                  wF => wC(:,:,this%kst:this%ken), fxF => this%fx(:,:,this%kst:this%ken), &
                  fyF => this%fy(:,:,this%kst:this%ken), fzF => rbuffxC(:,:,this%kst:this%ken,1), &
                  dudx => duidxjC(:,:,this%kst:this%ken,1), dudy => duidxjC(:,:,this%kst:this%ken,2), &
                  dudz => duidxjC(:,:,this%kst:this%ken,3), dvdx => duidxjC(:,:,this%kst:this%ken,4), &
                  dvdy => duidxjC(:,:,this%kst:this%ken,5), dvdz => duidxjC(:,:,this%kst:this%ken,6), &
                  dwdx => duidxjC(:,:,this%kst:this%ken,7), dwdy => duidxjC(:,:,this%kst:this%ken,8), &
                  dwdz => duidxjC(:,:,this%kst:this%ken,9), nSGS_F => nSGS(:,:,this%kst:this%ken))
              
                  call this%updateForce(rbuffxC, rbuffxE, cbuffxC, cbuffxE, cbuffyC, cbuffyE, cbuffzC, cbuffzE)

                  ! Step 3b: Time integrate (since dfdt = <Bodart forcing>) using Forwared Euler
                  this%fxhat = this%fxhat_old + dt*this%fxhat
                  this%fyhat = this%fyhat_old + dt*this%fyhat
                  this%fzhat = this%fzhat_old + dt*this%fzhat

                  ! Update and store fhat_old
                  this%fxhat_old = this%fxhat
                  this%fyhat_old = this%fyhat
                  this%fzhat_old = this%fzhat

                  ! Compare the current work done by the force to the target value
                  ! and scale the force appropriately
                  call this%spectC%ifft(this%fxhat,this%fx)
                  call this%spectC%ifft(this%fyhat,this%fy)
                  call this%spectE%ifft(this%fzhat,this%fz)

                  call this%interpE2C(this%fz,rbuffxC(:,:,:,1),rbuffyC(:,:,:,1),&
                    rbuffzC(:,:,:,1),rbuffyE(:,:,:,1),rbuffzE(:,:,:,1))
                  !forceWork    = p_sum(sum(u(:,:,this%kst:this%ken)*this%fx(:,:,this%kst:this%ken) + &
                  !  v(:,:,this%kst:this%ken)*this%fy(:,:,this%kst:this%ken) + &
                  !  wC(:,:,this%kst:this%ken)*rbuffxC(:,:,this%kst:this%ken,1)))

                  ! Now remove the mean if useing version 1 since version 2 is
                  ! guaranteed to have zero mean
                  if (this%spectC%carryingZeroK .and. (this%version == 1)) then
                      this%fxhat(this%spectC%zeroK_i,this%spectC%zeroK_j,:) = im0
                      this%fyhat(this%spectC%zeroK_i,this%spectC%zeroK_j,:) = im0
                  end if

                  forceWork    = p_sum(sum(uF*fxF + vF*fyF + wF*fzF))
                  select case (this%controlMethod)
                  case (1)
                    this%ampFact = this%tgtDissipation/(forceWork + 1.d-14)
                  case (2)
                    KE = 0.5d0*p_sum(sum(uF*uF + vF*vF + wF*wF))
                    eps = p_sum(sum((1.d0/Re + nSGS_F)*(dudx*dudx + dudy*dudy + &
                      dudz*dudz + dvdx*dvdx + dvdy*dvdy + dvdz*dvdz + &
                      dwdx*dwdx + dwdy*dwdy + dwdz*dwdz)))

                    this%ampFact = (eps - this%gain*(KE - this%tgtKE)/this%integralTime)/(forceWork + 1.d-14)
                  case (3) ! A(t) = const
                    this%ampFact = 1.d-2
                  end select

                  ! Clip large values due to forceWork~0
                  if (abs(this%ampFact) > 1.d4) this%ampFact = 1.d0
              end associate
          end if

          ! Step 4: Update RHS velocities
          urhs = urhs + this%ampFact*this%fxhat
          vrhs = vrhs + this%ampFact*this%fyhat
          wrhs = wrhs + this%ampFact*this%fzhat

          if (newTimeStep .and. this%checkDivergence) then
              call this%divergenceCheck(rbuffxC(:,:,:,1))
          end if

      end subroutine

      subroutine divergenceCheck(this,rbuffxC,fixDivergence)
          class(forcingLayer), intent(inout) :: this
          real(rkind), dimension(:,:,:), intent(inout) :: rbuffxC
          logical, intent(in), optional :: fixDivergence
          logical :: fixDivrg

          if (present(fixDivergence)) then
              fixDivrg = fixDivergence
          else
              fixDivrg = .false.
          end if

          rbuffxC = 0.d0
          call this%poiss%divergenceCheck(this%fxhat, this%fyhat, this%fzhat, &
            rbuffxC, fixDiv = fixDivrg, printMessage = .false.)
          this%maxDiv = p_maxval(maxval(abs(rbuffxC)))
          call message(1, "Max divergence of forcing layer",this%maxDiv)
      end subroutine
            
      subroutine updateForce(this, rbuffxC, rbuffxE, cbuffxC, cbuffxE, cbuffyC, cbuffyE, cbuffzC, cbuffzE)
          use random, only: uniform_random
          class(forcingLayer), intent(inout) :: this
          real(rkind), dimension(:,:,:,:), intent(inout) :: rbuffxC, rbuffxE
          complex(rkind), dimension(:,:,:,:), intent(inout) :: cbuffyC, cbuffyE, cbuffzC, cbuffzE
          complex(rkind), dimension(:,:,:), intent(inout) :: cbuffxC, cbuffxE
          integer :: i, j, n
          real(rkind), dimension(size(this%xs)) :: xmid, ymid
          real(rkind) :: xst, yst
          real(rkind) :: tmp
          real(rkind), dimension(this%nblocks,this%nblocks) :: sgn
          real(rkind), dimension(0:1) :: sgns
             
          this%fx = 0.d0
          this%fy = 0.d0
          this%fz = 0.d0
          
          this%xshift = 0.d0
          this%yshift = 0.d0
          this%zshift = 0.d0
          
          sgns = [-1.d0,1.d0]
          
          ! Step 1: get new locations for forcing blocks
          if (this%randomizeBlockPosition_xy) then
              call this%updateSeed()
              call uniform_random(this%xshift,-0.5d0*this%Lx,0.5d0*this%Lx,this%seed)
              call this%updateSeed()
              call uniform_random(this%yshift,-0.5d0*this%Ly,0.5d0*this%Ly,this%seed)
          end if
          if (this%randomizeBlockPosition_z) then
              ! Randomly shift each block in z-direction
              do j = 1,this%nblocks
                  do i = 1,this%nblocks
                      call this%updateSeed()
                      call uniform_random(this%zshift(i,j),-0.5d0*this%zgap,0.5d0*this%zgap,this%seed)
                  end do
              end do
          end if

          ! Elementary force types
          if (this%fixedForceType > 0) then
              this%forceType = this%fixedForceType
          else
              do j = 1,this%nblocks
                  do i = 1,this%nblocks
                      call this%updateSeed()
                      call uniform_random(tmp,0.d0,1.d0,this%seed)
                      this%forceType(i,j) = ceiling(tmp/third)

                      call this%updateSeed()
                      call uniform_random(tmp,0.d0,1.d0,this%seed)
                      sgn(i,j) = sgns(nint(tmp))
                  end do
              end do
          end if
          
          ! Step 2: get the elementary force for each block
          do j = 1,this%nblocks
              yst = (j - 1)*this%lf
              do i = 1,this%nblocks
                xst = (i - 1)*this%lf
                
                this%phixC = 0.d0; this%phiyC = 0.d0; this%phizC = 0.d0
                this%phixE = 0.d0; this%phiyE = 0.d0; this%phizE = 0.d0
                
                ! Add periodic contributions
                do n = 1,size(this%periodicXshift)
                    xmid = this%xs + xst + this%xshift + this%periodicXshift(n)
                    ymid = this%ys + yst + this%yshift + this%periodicYshift(n)
                
                    this%phixC = this%phixC + spline2(xmid,this%vs,this%x ,1,this%gpC)
                    this%phiyC = this%phiyC + spline2(ymid,this%vs,this%y ,2,this%gpC)
                    
                    this%phixE = this%phixE + spline2(xmid,this%vs,this%xE,1,this%gpE)
                    this%phiyE = this%phiyE + spline2(ymid,this%vs,this%yE,2,this%gpE)
                end do
                
                this%phizC = spline2(this%zs + this%zst + this%zshift(i,j),this%vs,this%z ,3,this%gpC)
                this%phizE = spline2(this%zs + this%zst + this%zshift(i,j),this%vs,this%zE,3,this%gpE)
               
                ! Get derivatives
                call ddx(this%phixC,this%dphixC,cbuffyC,this%spectC,this%sp_buffyC,this%k1C)
                call ddx(this%phixE,this%dphixE,cbuffyE,this%spectE,this%sp_buffyE,this%k1E)
                
                call ddy(this%phiyC,this%dphiyC,cbuffyC,this%spectC)
                call ddy(this%phiyE,this%dphiyE,cbuffyE,this%spectE)
                
                call ddz_E2C(this%phizE,this%dphizC,cbuffyE(:,:,:,1),cbuffyC(:,:,:,1),&
                  cbuffzE(:,:,:,1),cbuffzC(:,:,:,1), this%spectC,this%spectE,this%sp_gpE,&
                  this%sp_gpC,this%pade6opZ)
                call ddz_C2E(this%phizC,this%dphizE,cbuffyE(:,:,:,1),cbuffyC(:,:,:,1),&
                  cbuffzE(:,:,:,1),cbuffzC(:,:,:,1), this%spectC,this%spectE,this%sp_gpE,&
                  this%sp_gpC,this%pade6opZ)

                ! Randomly choose the force type
                call getForceBlock(this%phixC,this%phiyC,this%phizC, &
                  this%dphixC, this%dphiyC, this%dphizC, this%forceType(i,j), &
                  rbuffxC(:,:,:,1), rbuffxC(:,:,:,2),rbuffxC(:,:,:,3),this%version)
                call getForceBlock(this%phixE,this%phiyE,this%phizE, &
                  this%dphixE, this%dphiyE, this%dphizE, this%forceType(i,j), &
                  rbuffxE(:,:,:,1), rbuffxE(:,:,:,2), rbuffxE(:,:,:,3),this%version)

                this%fx = this%fx + sgn(i,j)*rbuffxC(:,:,:,1)
                this%fy = this%fy + sgn(i,j)*rbuffxC(:,:,:,2)
                this%fz = this%fz + sgn(i,j)*rbuffxE(:,:,:,3)

              end do
          end do
          ! Step 3: Take FFT (in x & y) of the force
          call this%spectC%fft(this%fx,this%fxhat)
          call this%spectC%fft(this%fy,this%fyhat)
          call this%spectE%fft(this%fz,this%fzhat)

          call fixOddBall(this%fxhat,this%spectC,cbuffxC)
          call fixOddBall(this%fyhat,this%spectC,cbuffxC)
          call fixOddBall(this%fzhat,this%spectE,cbuffxE)

      end subroutine

      subroutine ddz_E2C(f,dfdz,cbuffyE,cbuffyC,cbuffzE,cbuffzC,spectC,spectE,gpE,gpC,pade6opZ)
          real(rkind), dimension(:,:,:), intent(in) :: f
          real(rkind), dimension(:,:,:), intent(out) :: dfdz
          complex(rkind), dimension(:,:,:), intent(inout) :: cbuffyE, cbuffyC
          complex(rkind), dimension(:,:,:), intent(inout) :: cbuffzE, cbuffzC
          class(spectral), intent(inout) :: spectC, spectE
          class(decomp_info), intent(in) :: gpE, gpC
          class(Pade6stagg), intent(inout) :: pade6opZ

          call spectE%fft(f,cbuffyE)
          call transpose_y_to_z(cbuffyE,cbuffzE,gpE)
          call pade6opZ%ddz_E2C(cbuffzE, cbuffzC,zbc_bot,zbc_top)
          call transpose_z_to_y(cbuffzC,cbuffyC,gpC)
          call spectC%ifft(cbuffyC,dfdz)
      end subroutine
      
      subroutine ddz_C2E(f,dfdz,cbuffyE,cbuffyC,cbuffzE,cbuffzC,spectC,spectE,gpE,gpC,pade6opZ)
          real(rkind), dimension(:,:,:), intent(in) :: f
          real(rkind), dimension(:,:,:), intent(out) :: dfdz
          complex(rkind), dimension(:,:,:), intent(inout) :: cbuffyE, cbuffyC
          complex(rkind), dimension(:,:,:), intent(inout) :: cbuffzE, cbuffzC
          class(spectral), intent(inout) :: spectC, spectE
          class(decomp_info), intent(in) :: gpE, gpC
          class(Pade6stagg), intent(inout) :: pade6opZ

          call spectC%fft(f,cbuffyC)
          call transpose_y_to_z(cbuffyC,cbuffzC,gpC)
          call pade6opZ%ddz_C2E(cbuffzC,cbuffzE,zbc_bot,zbc_top)
          call transpose_z_to_y(cbuffzE,cbuffyE,gpE)
          call spectE%ifft(cbuffyE,dfdz)
      end subroutine

      subroutine ddx(f,dfdx,cbuffy,spect,rbuffy,myk1)
          real(rkind), dimension(:,:,:), intent(in) :: f
          real(rkind), dimension(:,:,:), intent(out) :: dfdx
          complex(rkind), dimension(:,:,:,:), intent(inout) :: cbuffy
          class(spectral), intent(inout) :: spect
          real(rkind), dimension(:,:,:), intent(inout) :: rbuffy
          real(rkind), dimension(:,:,:), intent(in) :: myk1

          call spect%fft(f,cbuffy(:,:,:,1))
          ! Zero the oddball wavenumber regardless if fixOdball=.false. in spect class
          rbuffy = spect%k1
          spect%k1 = myk1
          call spect%mTimes_ik1_oop(cbuffy(:,:,:,1),cbuffy(:,:,:,2))
          call spect%ifft(cbuffy(:,:,:,2), dfdx)
          spect%k1 = rbuffy
      end subroutine 
                
      subroutine ddy(f,dfdy,cbuffy,spect)
          real(rkind), dimension(:,:,:), intent(in) :: f
          real(rkind), dimension(:,:,:), intent(out) :: dfdy
          complex(rkind), dimension(:,:,:,:), intent(inout) :: cbuffy
          class(spectral), intent(inout) :: spect

          call spect%fft(f,cbuffy(:,:,:,1))
          cbuffy(:,spect%ny_g/2+1,:,1) = 0
          call spect%mTimes_ik2_oop(cbuffy(:,:,:,1),cbuffy(:,:,:,2))
          call spect%ifft(cbuffy(:,:,:,2), dfdy)
      end subroutine 

      subroutine fixOddball(fhat,spect,cbuffx)
          complex(rkind), dimension(:,:,:), intent(inout) :: fhat, cbuffx
          class(spectral), intent(inout) :: spect
          
          fhat(:,spect%ny_g/2+1,:) = 0
          call transpose_y_to_x(fhat,cbuffx,spect%spectdecomp)
          cbuffx(spect%nx_g/2+1,:,:) = 0
          call transpose_x_to_y(cbuffx,fhat,spect%spectdecomp)
      end subroutine
      
      subroutine interpC2E(this,fC,fE,rbuffyC,rbuffzC,rbuffyE,rbuffzE)
          class(forcingLayer), intent(inout) :: this
          real(rkind), dimension(:,:,:), intent(in) :: fC
          real(rkind), dimension(:,:,:), intent(out) :: fE
          real(rkind), dimension(:,:,:), intent(inout) :: rbuffyC, rbuffyE
          real(rkind), dimension(:,:,:), intent(inout) :: rbuffzC, rbuffzE
          call transpose_x_to_y(fC,rbuffyC,this%gpC)
          call transpose_y_to_z(rbuffyC, rbuffzC, this%gpC)
          call this%pade6opZ%interpz_C2E(rbuffzC, rbuffzE, zbc_bot, zbc_top)
          call transpose_z_to_y(rbuffzE, rbuffyE, this%gpE)
          call transpose_y_to_x(rbuffyE, fE, this%gpE)
      end subroutine
      
      subroutine interpE2C(this,fE,fC,rbuffyC,rbuffzC,rbuffyE,rbuffzE)
          class(forcingLayer), intent(inout) :: this
          real(rkind), dimension(:,:,:), intent(in) :: fE
          real(rkind), dimension(:,:,:), intent(out) :: fC
          real(rkind), dimension(:,:,:), intent(inout) :: rbuffyC, rbuffyE
          real(rkind), dimension(:,:,:), intent(inout) :: rbuffzC, rbuffzE
          
          call transpose_x_to_y(fE,rbuffyE,this%gpE)
          call transpose_y_to_z(rbuffyE,rbuffzE,this%gpE)
          call this%pade6opZ%interpz_E2C(rbuffzE, rbuffzC, zbc_bot, zbc_top)
          call transpose_z_to_y(rbuffzC, rbuffyC, this%gpC)
          call transpose_y_to_x(rbuffyC, fC, this%gpC)
      end subroutine
      
      subroutine updateSeed(this)
          use seedGen, only: incrementLogisticMap
          class(forcingLayer), intent(inout) :: this
          call incrementLogisticMap(this%seedFact,this%seedFact)
          this%seed = nint(rmaxInt*this%seedFact)
      end subroutine

      subroutine getForceBlock(phix,phiy,phiz,dphix,dphiy,dphiz,forceType,fx,fy,fz,version)
          real(rkind), dimension(:,:,:), intent(in) :: phix, phiy, phiz, dphix, &
            dphiy, dphiz
          integer, intent(in) :: forcetype
          real(rkind), dimension(:,:,:), intent(out) :: fx, fy, fz
          integer, intent(in) :: version
          
          fx = 0.d0; fy = 0.d0; fz = 0.d0

          select case (version)
          case (1)
              select case (forceType)
              case (1)
                  fx = phix*dphiy*phiz
                  fy = -dphix*phiy*phiz
                  fz = 0.d0
              case (2)
                  fx = -dphiz*phix*phiy
                  fy = 0.d0
                  fz = phiz*dphix*phiy
              case (3)
                  fx = 0.d0
                  fy = phiy*dphiz*phix
                  fz = -dphiy*phiz*phix
              end select
          case (2)
              select case (forceType)
              case (1)
                  fx = phix*dphiy*phiz*dphiz
                  fy = -dphix*phiy*phiz*dphiz
                  fz = 0.d0
              case (2)
                  fx = -dphiz*phix*phiy*dphiy
                  fy = 0.d0
                  fz = phiz*dphix*phiy*dphiy
              case (3)
                  fx = 0.d0
                  fy = phiy*dphiz*phix*dphix
                  fz = -dphiy*phiz*phix*dphix
              end select
          end select
      end subroutine
      
      function spline2(xs,vs,xq,coord,gp) result(vq)
          ! Second order spline interpolation
          real(rkind), dimension(:), intent(in) :: xs, vs
          real(rkind), dimension(:,:,:), intent(in) :: xq
          integer, intent(in) :: coord
          class(decomp_info), intent(in) :: gp
          real(rkind), dimension(gp%xsz(1),gp%xsz(2),gp%xsz(3)) :: vq
          real(rkind), dimension(size(xq,coord)) :: x1D
          real(rkind), dimension(size(xs)) :: weight
          integer :: nknots, st, en, n
          real(rkind) :: delta

          select case (coord)
          case (1)
            x1D = xq(:,1,1)
          case (2)
            x1D = xq(1,:,1)
          case (3)
            x1D = xq(1,1,:)
          end select

          ! The number of sample points (i.e. "knots")
          nknots = size(xs)

          ! Assumes uniform grid spacing
          delta = xs(2) - xs(1)

          ! Compute the weights (g-primes)
          weight = getWeights(vs,0.d0,delta)

          vq = 0.d0
          
          ! Now assemble the spline function
          select case (coord)
          case (1)
            do n = 1,nknots-1
              if (onThisRank(x1D,xs(n),xs(n+1))) then
                  call getStEndIndices(x1D,xs(n),xs(n+1),st,en)
                  vq(st:en,:,:) = spl2fcn(xq(st:en,:,:),vs(n),weight(n),weight(n+1),delta,xs(n),xs(n+1))
              end if
            end do
          case (2)
            do n = 1,nknots-1
              if (onThisRank(x1D,xs(n),xs(n+1))) then
                  call getStEndIndices(x1D,xs(n),xs(n+1),st,en)
                  vq(:,st:en,:) = spl2fcn(xq(:,st:en,:),vs(n),weight(n),weight(n+1),delta,xs(n),xs(n+1))
              end if
            end do
          case (3)
            do n = 1,nknots-1
              if (onThisRank(x1D,xs(n),xs(n+1))) then
                  call getStEndIndices(x1D,xs(n),xs(n+1),st,en)
                  vq(:,:,st:en) = spl2fcn(xq(:,:,st:en),vs(n),weight(n),weight(n+1),delta,xs(n),xs(n+1))
              end if
            end do
          end select
      end

      function spl2fcn(xq,f1,gp1,gp2,delta,x1,x2) result (vq)
          real(rkind), dimension(:,:,:), intent(in) :: xq
          real(rkind), intent(in) :: f1, gp1, gp2, delta, x1, x2
          real(rkind), dimension(size(xq,1),size(xq,2),size(xq,3)) :: vq
          real(rkind) :: C

          C = f1 + delta/2.d0*gp1
          vq = -gp1/(2.d0*delta)*(xq - x2)*(xq - x2) + gp2/(2.d0*delta)*(xq - x1)*(xq - x1) + C

      end function

      function getWeights(f,gp0,delta) result (gp)
          real(rkind), dimension(:), intent(in) :: f
          real(rkind), intent(in) :: gp0, delta
          real(rkind), dimension(size(f)) :: gp
          integer :: i, n

          n = size(f)
          gp(1) = gp0
          do i = 1,n-1
            gp(i+1) = 2*(f(i+1) - f(i))/delta - gp(i)
          end do
      end function

      function onThisRank(myx,xst,xen) result (TF)
          real(rkind), dimension(:), intent(in) :: myx
          real(rkind), intent(in) :: xst, xen
          logical :: TF

          TF = .true.
          if ((xen < minval(myx)) .or. (xst > maxval(myx))) TF = .false.
      end function
      
      subroutine getStEndIndices(xvec,x1,x2,st,en)
          real(rkind), dimension(:), intent(in) :: xvec
          real(rkind), intent(in) :: x1, x2
          integer, intent(out) :: st, en

          st = minloc(abs(xvec - x1), dim=1)
          en = minloc(abs(xvec - x2), dim=1)
          if (xvec(st) < x1) st = st + 1
          if (xvec(en) > x2) en = en - 1
      end subroutine
   
end module
