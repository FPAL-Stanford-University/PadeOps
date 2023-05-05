module forcingLayerMod
  ! This forcing scheme is adapted from that presented in 
  ! Bodart, Julien, J-B. Cazalbou, and Laurent Joly. "Direct numerical
  ! simulation of unsheared turbulence diffusing towards a free-slip or no-slip
  ! surface." Journal of Turbulence 11 (2010): N48.
    use kind_parameters, only: rkind, clen, mpirkind
    use decomp_2d, only: decomp_info, transpose_x_to_y, transpose_y_to_z, &
        transpose_z_to_y, transpose_y_to_x, nrank
    use spectralMod, only: spectral
    use reductions, only: p_minval, p_maxval
    use PadeDerOps, only: Pade6Stagg
    use seedGen, only: initializeLogMap
    use PadePoissonMod, only: Padepoisson 
    use exits, only: message
    use basic_io, only: read_0d_ascii
    use mpi
    
    implicit none
    private

    real(rkind), dimension(4) :: xsample = [0.d0, 1.d0/3.d0, 2.d0/3.d0, 1.d0]
    real(rkind), dimension(4) :: vsample = [0.d0, 2.d0/3.d0, 2.d0/3.d0, 0.d0]
    integer, parameter :: zbc_bot = -1, zbc_top = -1
    public :: forcingLayer

    type :: forcingLayer
        ! TODO: 
        !   [x] Dump seed as restart file (seed is the only information needed)
        !   [_] Implement periodicity 
        integer :: nblocks ! The number of forcing blocks in each x,y direction
        real(rkind) :: lf  ! Physical size of forcing layer
        real(rkind), dimension(:,:,:), allocatable :: fx, fy, fz
        real(rkind), dimension(:,:,:), allocatable :: phixC, phiyC, phizC
        real(rkind), dimension(:,:,:), allocatable :: phixE, phiyE, phizE
        real(rkind), dimension(:,:,:), allocatable :: dphixC, dphiyC, dphizC
        real(rkind), dimension(:,:,:), allocatable :: dphixE, dphiyE, dphizE
        type(decomp_info), pointer :: gpC, gpE, sp_gpC, sp_gpE
        type(spectral), pointer :: spectC, spectE
        real(rkind), dimension(4) :: xs, ys, zs, vs ! Sample points for spline interpolation
        real(rkind), dimension(:,:,:), pointer :: x, y, z, zE
        real(rkind) :: zgap, forceAmp
        type(Pade6Stagg), pointer :: Pade6opZ
        logical :: dumpForce, dumpSplines, checkDivergence
        logical :: randomizeBlockPositions
        real(rkind) :: seedFact
        integer :: seed, fixedForceType
        real(rkind) :: xshift = 0.d0, yshift = 0.d0
        real(rkind), dimension(:,:), allocatable :: zshift
        type(padePoisson), pointer :: poiss
        integer, dimension(:,:), allocatable :: forceType
        real(rkind) :: Lx, Ly
        real(rkind), dimension(:,:,:), allocatable :: k1C, k1E
        real(rkind), dimension(:,:,:), allocatable :: sp_buffyC, sp_buffyE
        contains
          procedure :: init
          procedure :: destroy
          procedure :: updateRHS
          procedure :: interpC2E
          procedure :: interpE2C
          procedure :: updateSeed
    end type

    contains

      subroutine init(this,inputfile,Lx,Ly,mesh,zE,gpC,gpE,spectC,spectE,Pade6opZ,PadePoiss,&
          restartSim,runID,tsim,inputdir)
          class(forcingLayer), intent(inout) :: this
          character(len=*), intent(in) :: inputfile, inputdir
          real(rkind), intent(in) :: Lx, Ly
          real(rkind), dimension(:,:,:,:), intent(in), target :: mesh ! Cell-centered mesh (x,y,z)
          real(rkind), dimension(:,:,:), intent(in), target :: zE ! Edge mesh (z only)
          class(decomp_info), intent(in), target :: gpC, gpE
          class(spectral), intent(in), target :: spectC, spectE
          class(Pade6Stagg), intent(in), target :: Pade6opZ
          class(padePoisson), intent(in), target :: PadePoiss
          logical, intent(in) :: restartSim
          integer, intent(in) :: runID, tsim
          integer :: nblocks = 1, i, j
          integer :: ioUnit, ierr, fixedForceType = 0
          real(rkind) :: zgap = 0, forceAmp = 1
          logical :: dumpForce = .false., dumpSplines = .false., checkDivergence = .false.
          logical :: randomizeBlockPositions = .true.
          character(len=clen) :: tempname, fname
          real(rkind), dimension(:,:,:), allocatable :: tmp

          namelist /localizedForceLayer/ nblocks, zgap, forceAmp, dumpForce, &
            dumpSplines, checkDivergence, fixedForceType, randomizeBlockPositions
          
          ioUnit = 123
          open(unit=ioUnit, file=trim(inputfile), form='FORMATTED', iostat=ierr)
          read(unit=ioUnit, NML=localizedForceLayer)
          close(ioUnit)

          ! Input file parameters
          !     nblocks --> Number of forcing blocks in each x, y direction
          !     forceAmp --> Amplitude of the forcing
          !               we can randomly shift the forcing blocks in x and y without worrying
          !               about periodicity
          !     zgap --> Allows random shifts in z

          this%nblocks  = nblocks 
          this%forceAmp = forceAmp 
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
          this%vs = vsample!*this%forceAmp
          
          this%x => mesh(:,:,:,1)
          this%y => mesh(:,:,:,2)
          this%z => mesh(:,:,:,3)
          this%zE => zE

          this%dumpForce = dumpForce
          this%dumpSplines = dumpSplines
          this%checkDivergence = checkDivergence
          this%randomizeBlockPositions = randomizeBlockPositions

          this%fixedForceType = fixedForceType
        
          allocate(this%forceType(this%nblocks,this%nblocks))

          allocate(this%fx(this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3)))
          allocate(this%fy(this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3)))
          allocate(this%fz(this%gpE%xsz(1), this%gpE%xsz(2), this%gpE%xsz(3)))
          
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

          ! Initialize randomized seeds
          if (restartSim) then
              if (nrank == 0) then
                  write(tempname,"(A7,A4,I2.2,A5,I6.6)") "RESTART", "_Run",runID, "_frc.",tsim
                  fname = inputDir(:len_trim(inputDir))//"/"//trim(tempname)
                  call read_0d_ascii(this%seedFact,fname)
              end if
              call MPI_Barrier(MPI_COMM_WORLD,ierr)
              call MPI_Bcast(this%seedFact,1,mpirkind,0,MPI_COMM_WORLD,ierr)
              call MPI_Barrier(MPI_COMM_WORLD,ierr)
          else
              this%seedFact = initializeLogMap(1)
          end if

          call message(0,'Forcing layer initialized')
          call message(1,'Forcing layer bounds (about forcing mid-plane)',[-0.5d0*(this%zgap + this%lf), &
            0.5d0*(this%zgap + this%lf)])
      end subroutine

      subroutine destroy(this)
          class(forcingLayer), intent(inout) :: this

          deallocate(this%fx,this%fy,this%fz)
          deallocate(this%phixC,this%phiyC,this%phizC)
          deallocate(this%phixE,this%phiyE,this%phizE)
          deallocate(this%dphixC,this%dphiyC,this%dphizC)
          deallocate(this%dphixE,this%dphiyE,this%dphizE)
          deallocate(this%zshift)
          deallocate(this%forceType)

          nullify(this%gpC, this%gpE, this%spectC, this%spectE, this%sp_gpC)
          nullify(this%sp_gpE, this%Pade6opz, this%x, this%y, this%z, this%poiss)
          nullify(this%zE)
      end subroutine

      subroutine updateRHS(this,urhs,vrhs,wrhs,cbuffxC,cbuffxE,cbuffyC,cbuffyE,cbuffzC,cbuffzE,&
          rbuffxC,rbuffxE,rbuffyC,rbuffzC,rbuffyE,rbuffzE,newTimeStep)
          use random, only: uniform_random
          use constants, only: third, half
          class(forcingLayer), intent(inout) :: this
          complex(rkind), dimension(this%sp_gpC%ysz(1), this%sp_gpC%ysz(2), this%sp_gpC%ysz(3)), intent(inout) :: urhs, vrhs
          complex(rkind), dimension(this%sp_gpE%ysz(1), this%sp_gpE%ysz(2), this%sp_gpE%ysz(3)), intent(inout) :: wrhs
          complex(rkind), dimension(:,:,:), intent(inout) :: cbuffxC, cbuffxE
          complex(rkind), dimension(:,:,:,:), intent(inout) :: cbuffyC, cbuffyE
          real(rkind), dimension(:,:,:), intent(inout) :: rbuffyC, rbuffzC, rbuffyE, rbuffzE
          real(rkind), dimension(:,:,:,:), intent(inout) :: rbuffxC, rbuffxE
          logical, intent(in) :: newTimeStep
          real(rkind) :: xst, yst, zst
          integer :: i, j, n
          real(rkind) :: tmp, ampFact = 1.d0
          complex(rkind), dimension(:,:,:,:), intent(inout) :: cbuffzC, cbuffzE
          real(rkind), dimension(size(this%xs)) :: xmid, ymid
          real(rkind), dimension(3) :: periodicXshift, periodicYshift
          real(rkind), dimension(this%nblocks,this%nblocks) :: sgn
          real(rkind), dimension(0:1) :: sgns

          sgns = [-1.d0,1.d0]

          ! Step 1: get new locations for forcing blocks
          if (newTimeStep) then
              if (this%randomizeBlockPositions) then
                  call this%updateSeed()
                  call uniform_random(this%xshift,-0.5d0*this%Lx,0.5d0*this%Lx,this%seed)
                  call this%updateSeed()
                  call uniform_random(this%yshift,-0.5d0*this%Ly,0.5d0*this%Ly,this%seed)
                  ! Randomly shift each block in z-direction
                  do j = 1,this%nblocks
                      do i = 1,this%nblocks
                          call this%updateSeed()
                          call uniform_random(this%zshift(i,j),-0.5d0*this%zgap,0.5d0*this%zgap,this%seed)
                      end do
                  end do
              else
                  this%xshift = 0.d0
                  this%yshift = 0.d0
                  this%zshift = 0.d0
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
          end if
          ! Location of the bottom of the forcing layer
          zst = -0.5d0*this%lf

          periodicXshift = [0.d0, this%Lx, -this%Lx]
          periodicYshift = [0.d0, this%Ly, -this%Ly]

          if (newTimeStep) then
              this%fx = 0.d0
              this%fy = 0.d0
              this%fz = 0.d0
              
              ! Step 2: get the elementary force for each block
              do j = 1,this%nblocks
                  yst = (j - 1)*this%lf
                  do i = 1,this%nblocks
                    xst = (i - 1)*this%lf
                    
                    this%phixC = 0.d0; this%phiyC = 0.d0; this%phizC = 0.d0
                    this%phixE = 0.d0; this%phiyE = 0.d0; this%phizE = 0.d0
                    
                    ! Loop over all 8 periodic shifts
                    do n = 1,size(periodicXshift)
                        xmid = this%xs + xst + this%xshift + periodicXshift(n)
                        ymid = this%ys + yst + this%yshift + periodicYshift(n)
                    
                        this%phixC = this%phixC + spline2(xmid,this%vs,this%x ,1,this%gpC)
                        this%phiyC = this%phiyC + spline2(ymid,this%vs,this%y ,2,this%gpC)
                        
                        this%phixE = this%phixE + spline2(xmid,this%vs,this%x ,1,this%gpE)
                        this%phiyE = this%phiyE + spline2(ymid,this%vs,this%y ,2,this%gpE)
                    end do
                    
                    this%phizC = spline2(this%zs + zst + this%zshift(i,j),this%vs,this%z ,3,this%gpC)
                    this%phizE = spline2(this%zs + zst + this%zshift(i,j),this%vs,this%zE,3,this%gpE)
                   
                    ! Get derivatives
                    call ddx(this%phixC,this%dphixC,cbuffyC,cbuffxC,this%spectC,this%sp_buffyC,this%k1C)
                    call ddx(this%phixE,this%dphixE,cbuffyE,cbuffxE,this%spectE,this%sp_buffyE,this%k1E)
                    
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
                      rbuffxC(:,:,:,1), rbuffxC(:,:,:,2), rbuffxC(:,:,:,3))
                    call getForceBlock(this%phixE,this%phiyE,this%phizE, &
                      this%dphixE, this%dphiyE, this%dphizE, this%forceType(i,j), &
                      rbuffxE(:,:,:,1), rbuffxE(:,:,:,2), rbuffxE(:,:,:,3))

                    !call this%updateSeed()
                    !call uniform_random(ampFact,0.d0,2.d0,this%seed)
                    this%fx = this%fx + ampFact*this%forceAmp*sgn(i,j)*rbuffxC(:,:,:,1)
                    this%fy = this%fy + ampFact*this%forceAmp*sgn(i,j)*rbuffxC(:,:,:,2)
                    this%fz = this%fz + ampFact*this%forceAmp*sgn(i,j)*rbuffxE(:,:,:,3)

                  end do
              end do
          end if

          ! Step 3: Take FFT (in x & y) of the force
          call this%spectC%fft(this%fx,cbuffyC(:,:,:,1))
          call this%spectC%fft(this%fy,cbuffyC(:,:,:,2))
          call this%spectE%fft(this%fz,cbuffyE(:,:,:,1))

          call fixOddBall(cbuffyC(:,:,:,1),cbuffyC(:,:,:,2),this%spectC,cbuffxC)

          ! Step 4: Update RHS velocities
          urhs = urhs + cbuffyC(:,:,:,1)
          vrhs = vrhs + cbuffyC(:,:,:,2)
          wrhs = wrhs + cbuffyE(:,:,:,1)

          if (newTimeStep .and. this%checkDivergence) then
              call this%poiss%divergenceCheck(cbuffyC(:,:,:,1),cbuffyC(:,:,:,2), &
                cbuffyE(:,:,:,1), rbuffxC(:,:,:,1))
              call message("Max divergence of forcing layer",p_maxval(maxval(abs(rbuffxC(:,:,:,1)))))
          end if

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

      subroutine ddx(f,dfdx,cbuffy,cbuffx,spect,rbuffy,myk1)
          real(rkind), dimension(:,:,:), intent(in) :: f
          real(rkind), dimension(:,:,:), intent(out) :: dfdx
          complex(rkind), dimension(:,:,:), intent(inout) :: cbuffx
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

      subroutine fixOddball(fxhat,fyhat,spect,cbuffx)
          complex(rkind), dimension(:,:,:), intent(inout) :: fxhat, fyhat, cbuffx
          class(spectral), intent(inout) :: spect
          
          fyhat(:,spect%ny_g/2+1,:) = 0
          call transpose_y_to_x(fxhat,cbuffx,spect%spectdecomp)
          cbuffx(spect%nx_g/2+1,:,:) = 0
          call transpose_x_to_y(cbuffx,fxhat,spect%spectdecomp)
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
          call transpose_y_to_z(rbuffyE, rbuffzE, this%gpE)
          call this%pade6opZ%interpz_E2C(rbuffzE, rbuffzC, zbc_bot, zbc_top)
          call transpose_z_to_y(rbuffzC, rbuffyC, this%gpC)
          call transpose_y_to_x(rbuffyC, fC, this%gpC)
      end subroutine
      
      subroutine updateSeed(this)
          use seedGen, only: incrementLogisticMap
          use constants, only: rmaxInt
          class(forcingLayer), intent(inout) :: this
          call incrementLogisticMap(this%seedFact,this%seedFact)
          this%seed = nint(rmaxInt*this%seedFact)
      end subroutine

      subroutine getForceBlock(phix,phiy,phiz,dphix,dphiy,dphiz,forceType,fx,fy,fz)
          real(rkind), dimension(:,:,:), intent(in) :: phix, phiy, phiz, dphix, &
            dphiy, dphiz
          integer, intent(in) :: forcetype
          real(rkind), dimension(:,:,:), intent(out) :: fx, fy, fz
          
          fx = 0.d0; fy = 0.d0; fz = 0.d0

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
          integer :: nknots, st, en
          real(rkind) :: delta, n

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
