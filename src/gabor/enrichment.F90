module enrichmentMod
  use kind_parameters,    only: rkind, clen, single_kind, castSingle, mpirkind
  use incompressibleGrid, only: igrid
  use QHmeshMod,          only: QHmesh
  use exits,              only: message
  use constants,          only: pi
  use fortran_assert,     only: assert
  use omp_lib
  use gaborHooks
  use timer,              only: tic, toc
  use decomp_2d,          only: DECOMP_2D_COMM_CART_X, nrank, nproc
  use mpi
  use reductions,         only: p_maxval, p_minval
  use timer,              only: tic, toc
  implicit none

  private
  integer :: nthreads
  integer, dimension(4) :: neighbour
  integer, dimension(2) :: coords, dims
  logical, dimension(3) :: periodicBCs

  real(rkind), dimension(2) :: xDom, yDom, zDom

  public :: enrichmentOperator, xDom, yDom, zDom, nthreads

  type :: enrichmentOperator
    !private
    real(rkind), dimension(:,:), allocatable :: modeData
    real(rkind), dimension(:), pointer, public :: kx, ky, kz
    real(rkind), dimension(:), pointer, public :: x, y, z
    real(rkind), dimension(:), pointer, public :: uhatR, uhatI, vhatR, &
      vhatI, whatR, whatI
    real(rkind) :: dt

    ! Halo-padded arrays
    real(rkind), dimension(:,:,:), allocatable :: uh, vh, wh
    real(rkind), dimension(:,:,:), allocatable :: KE, L, KEh, Lh
    real(rkind), dimension(:,:,:,:), allocatable :: duidxj_h

    integer, public :: nxsupp, nysupp, nzsupp
    type(igrid), pointer :: largeScales, smallScales
    type(QHmesh), public :: QHgrid 
    real(rkind), dimension(:,:,:,:), pointer :: duidxj_LS
    real(rkind) :: kmin, kmax
    logical :: imposeNoPenetrationBC = .false.

    ! Initialization parameters
    integer     :: nk, ntheta, nmodes, nmodesGlobal
    real(rkind) :: scalefact, Anu, numolec, ctauGlobal
    logical     :: renderPressure = .false.
    integer     :: tidRender, tio , tid, tidStop, tidsim

    ! Data IO
    character(len=clen) :: outputdir
    logical :: writeIsotropicModes
    
    ! Extra memory for velocity rendering
    real(single_kind), dimension(:,:,:,:), allocatable :: utmp,vtmp,wtmp
      ! Store info of modes on neighbor ranks that are within the support 
      ! window "halo"
    real(rkind), dimension(:,:), allocatable :: renderModeData

    ! Misc
    logical :: debugChecks = .false.
    logical :: strainInitialCondition = .true.

    contains
      procedure          :: init
      procedure          :: destroy
      procedure          :: generateIsotropicModesOld
      procedure          :: generateIsotropicModes
      procedure          :: strainModes
      procedure          :: getLargeScaleDataAtModeLocation
      procedure          :: advanceTime
      procedure          :: renderVelocity
      procedure          :: renderLocalVelocity
      procedure          :: updateLargeScales
      procedure          :: wrapupTimeStep
      procedure          :: dumpSmallScales
      procedure          :: continueSimulation
      procedure          :: dumpData
      procedure, private :: doDebugChecks
      procedure, private :: updateSeeds

      ! MPI communication stuff
      procedure          :: sendRecvHaloModes
  end type

contains
  include 'enrichment_files/initializeGaborModes.F90'
  include 'enrichment_files/renderVelocity.F90'
  include 'enrichment_files/gaborIO.F90'
  include 'enrichment_files/debugChecks.F90'
  include 'enrichment_files/MPIstuff.F90'

  subroutine init(this,smallScales,largeScales,inputfile)
    use GaborModeRoutines, only: computeKminKmax
    
    class(enrichmentOperator), intent(inout), target :: this 
    class(igrid), intent(inout), target :: smallScales, largeScales
    character(len=*), intent(in) :: inputfile
    character(len=clen) :: outputdir
    integer :: ierr, ioUnit
    integer :: nk, ntheta
    integer :: tidRender, tio, tidStop, tidInit = 0
    real(rkind) :: scalefact = 1.d0, Anu = 1.d-4, numolec = 0.d0
    real(rkind) :: ctauGlobal = 1.d0
    logical :: writeIsotropicModes = .false.
    integer :: ist, ien, jst, jen, kst, ken
    real(rkind) :: dt
    logical :: debugChecks = .false.
    logical :: strainInitialCondition = .true.
    logical :: xPeriodic = .true., yPeriodic = .true., zPeriodic = .true.
    integer :: n
    
    namelist /IO/      outputdir, writeIsotropicModes
    namelist /GABOR/   nk, ntheta, scalefact, ctauGlobal, Anu, numolec, &
      strainInitialCondition, xPeriodic, yPeriodic, zPeriodic
    namelist /CONTROL/ tidRender, tio, tidStop, tidInit, debugChecks
    namelist /INPUT/ dt

    ! Read inputfile
    ioUnit = 1
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED', iostat=ierr)
    read(unit=ioUnit, NML=INPUT)
    read(unit=ioUnit, NML=IO)
    read(unit=ioUnit, NML=GABOR)
    read(unit=ioUnit, NML=CONTROL)
    close(ioUnit)

    this%largeScales => largeScales 
    this%smallScales => smallScales

    ! Time control
    this%tidRender = tidRender 
    this%tio = tio
    this%tidStop = tidStop
    this%tid = tidInit
    this%dt = dt 
    this%debugChecks = debugChecks

    ! Initialization
    this%nk = nk
    this%ntheta = ntheta
    this%scalefact = scalefact
    this%Anu = Anu 
    this%ctauGlobal = ctauGlobal
    this%strainInitialCondition = strainInitialCondition

    ! IO
    this%outputdir = outputdir
    this%writeIsotropicModes = writeIsotropicModes

    ! Get domain boundaries
    xDom(1) = p_minval(this%largeScales%mesh(1,1,1,1))
    xDom(2) = p_maxval(this%largeScales%mesh(this%largeScales%gpC%xsz(1),1,1,1)) &
      + this%largeScales%dx
    yDom(1) = p_minval(this%largeScales%mesh(1,1,1,2))
    yDom(2) = p_maxval(this%largeScales%mesh(1,this%largeScales%gpC%xsz(2),1,2)) &
      + this%largeScales%dy
    zDom(1) = p_minval(this%largeScales%mesh(1,1,1,3)) &
      - 0.5d0*this%largeScales%dz
    zDom(2) = p_maxval(this%largeScales%mesh(1,1,this%largeScales%gpC%xsz(3),3)) &
      + 0.5d0*this%largeScales%dz

    ! STEP  0: Get filenames for u, v and w from (input file?) for initialization
    ! Use the same file-naming convention for PadeOps restart files: 
    ! RESTART_Run00_u.000000
    call this%largeScales%initLargeScales(tidInit,this%largeScales%runID)
    this%duidxj_LS => largeScales%duidxjC

    ! Get halo-padded velocity arrays
    allocate(this%KE(this%largeScales%gpC%xsz(1),this%largeScales%gpC%xsz(2),&
      this%largeScales%gpC%xsz(3)))
    allocate(this%L(this%largeScales%gpC%xsz(1),this%largeScales%gpC%xsz(2),&
      this%largeScales%gpC%xsz(3)))
    call getLargeScaleParams(this%KE,this%L,this%largeScales)
    call this%updateLargeScales(timeAdvance=.false.) 

    ! Ryan 
    ! STEP 1: Finish the QH region code to fill Gabor Modes (kx, ky, kz, x, y, z, uhat, ...)
    call this%QHgrid%init(inputfile,this%largeScales)

    this%nxsupp = 2 * this%smallScales%nx/this%largeScales%nx * &
      nint(this%QHgrid%dx / this%largeScales%dx)
    this%nysupp = 2 * this%smallScales%ny/this%largeScales%ny * &
      nint(this%QHgrid%dy / this%largeScales%dy)
    this%nzsupp = 2 * this%smallScales%nz/this%largeScales%nz * &
      nint(this%QHgrid%dz / this%largeScales%dz)

    ! User safegaurds -- make sure the mode supports don't span more than two
    ! MPI ranks
    call assert(p_minval(this%smallScales%gpC%xsz(2)) >= this%nysupp,&
      'The mode support widths span more than one MPI rank. Use fewer ranks'//&
      ' or adjust the topology -- xsz(2) >= nysupp')
    call assert(p_minval(this%smallScales%gpC%xsz(3)) >= this%nzsupp,&
      'The mode support widths span more than one MPI rank. Use fewer ranks'//&
      ' or adjust the topology -- xsz(3) >= nzsupp')

    ! Compute the number of modes
    this%nmodes = this%nk*this%ntheta * &
      this%QHgrid%gpC%xsz(1)*this%QHgrid%gpC%xsz(2)*this%QHgrid%gpC%xsz(3)
    this%nmodesGlobal = this%nk*this%ntheta * &
      this%QHgrid%nx*this%QHgrid%ny*this%QHgrid%nz
    
    ! Compute kmin and kmax based on LES and high-resolution grids 
    call computeKminKmax(xDom(2)-xDom(1), yDom(2)-yDom(1), zDom(2)-zDom(1), &
      this%largeScales%nx, this%largeScales%ny, this%largeScales%nz, &
      this%smallScales%nx, this%smallScales%ny, this%smallScales%nz, &
      this%kmin, this%kmax)

    ! Allocate memory for mode data
    allocate(this%modeData(this%nmodes,12))

    ! Link pointers
    this%x     => this%modeData(:,1)
    this%y     => this%modeData(:,2)
    this%z     => this%modeData(:,3)
    this%kx    => this%modeData(:,4)
    this%ky    => this%modeData(:,5)
    this%kz    => this%modeData(:,6)
    this%uhatR => this%modeData(:,7)
    this%uhatI => this%modeData(:,8)
    this%vhatR => this%modeData(:,9)
    this%vhatI => this%modeData(:,10)
    this%whatR => this%modeData(:,11)
    this%whatI => this%modeData(:,12)
   
    ! Allocate extra memory for velocity rendering
    ist = this%smallScales%gpC%xst(1) 
    ien = this%smallScales%gpC%xen(1) 
    jst = this%smallScales%gpC%xst(2) 
    jen = this%smallScales%gpC%xen(2) 
    kst = this%smallScales%gpC%xst(3) 
    ken = this%smallScales%gpC%xen(3)
   
    nthreads = omp_get_num_threads() 
    allocate(this%utmp(ist:ien,jst:jen,kst:ken,nthreads))
    allocate(this%vtmp(ist:ien,jst:jen,kst:ken,nthreads))
    allocate(this%wtmp(ist:ien,jst:jen,kst:ken,nthreads))

    ! Set things up for distributed memory
    call getNeighbours(neighbour)
    call MPI_Cart_Get(DECOMP_2D_COMM_CART_X,2,dims,periodicBCs(2:3),coords,ierr)
    !call this%largeScales%getPeriodicBCs(periodicBCs)
    periodicBCs(1) = xPeriodic
    periodicBCs(2) = yPeriodic
    periodicBCs(3) = zPeriodic

    ! Initialize the Gabor modes
    call this%generateIsotropicModesOld()
    !call this%generateIsotropicModes()

    if (this%writeIsotropicModes) then
      call message(1, 'Writing modes to disk.')
      call this%dumpData(this%x,this%y,this%z,this%kx,this%ky,this%kz, &
        this%uhatR,this%uhatI,this%vhatR,this%vhatI,this%whatR,this%whatI)
    end if
    if (this%strainInitialCondition) call this%strainModes()
   
    call this%wrapupTimeStep()
  end subroutine

  subroutine destroy(this)
    class(enrichmentOperator), intent(inout) :: this 

    if (associated(this%largeScales)) nullify(this%largeScales)
    if (associated(this%smallScales)) nullify(this%smallScales)
    if (allocated(this%modeData))    deallocate(this%modeData)
    if (associated(this%uhatR))    nullify(this%uhatR)
    if (associated(this%uhatI))    nullify(this%uhatI)
    if (associated(this%vhatR))    nullify(this%vhatR)
    if (associated(this%vhatI))    nullify(this%vhatI)
    if (associated(this%whatR))    nullify(this%whatR)
    if (associated(this%whatI))    nullify(this%whatI)
    if (associated(this%kx))       nullify(this%kx)
    if (associated(this%ky))       nullify(this%ky)
    if (associated(this%kz))       nullify(this%kz)
    if (associated(this%x))        nullify(this%x)
    if (associated(this%y))        nullify(this%y)
    if (associated(this%z))        nullify(this%z)
    if (allocated(this%uh))       deallocate(this%uh)
    if (allocated(this%vh))       deallocate(this%vh)
    if (allocated(this%wh))       deallocate(this%wh)
    if (allocated(this%duidxj_h)) deallocate(this%duidxj_h)
    if (allocated(this%utmp))     deallocate(this%utmp)
    if (allocated(this%vtmp))     deallocate(this%vtmp)
    if (allocated(this%wtmp))     deallocate(this%wtmp)

    call this%QHgrid%destroy()
  end subroutine 

  subroutine updateLargeScales(this,timeAdvance)
    class(enrichmentOperator), intent(inout) :: this 
    logical, intent(in) :: timeAdvance
    integer :: i

    ! Ryan: 
    ! STEP 1: Figure out how you want to advance the large-scales
    ! Either run a PadeOps time-step or, read in from a file 
    ! If you read in from a file, you would want to avoid doing this repeatedly. 
    ! Perhaps read in 20 snapshots at a time and so on..
    if (timeAdvance) call this%largeScales%timeAdvance(this%dt)

    call this%largeScales%fixGradientsForPeriodicity(periodicBCs)

    ! Aditya: 
    ! STEP 2: Generate halo'd velocities
    if(.not. allocated(this%uh))       call this%largeScales%pg%alloc_array(this%uh)
    if(.not. allocated(this%vh))       call this%largeScales%pg%alloc_array(this%vh)
    if(.not. allocated(this%wh))       call this%largeScales%pg%alloc_array(this%wh)
    if(.not. allocated(this%duidxj_h)) call this%largeScales%pg%alloc_array(this%duidxj_h,9)
    if(.not. allocated(this%KEh)) call this%largeScales%pg%alloc_array(this%KEh)
    if(.not. allocated(this%Lh)) call this%largeScales%pg%alloc_array(this%Lh)
    call this%largeScales%HaloUpdateVelocities(this%uh, this%vh, this%wh, &
      this%duidxj_h)

    call this%largeScales%haloUpdateField(this%KE,this%KEh)
    call this%largeScales%haloUpdateField(this%L, this%Lh)

    call removePeriodicityFromGhost(this%uh,this%largeScales%gpC)
    call removePeriodicityFromGhost(this%vh,this%largeScales%gpC)
    call removePeriodicityFromGhost(this%wh,this%largeScales%gpC)
    do i = 1,9
      call removePeriodicityFromGhost(this%duidxj_h(:,:,:,i),this%largeScales%gpC)
    end do
  end subroutine

  subroutine removePeriodicityFromGhost(f,gp)
    use decomp_2d, only: decomp_info
    use procgrid_mod, only: num_pad
    real(rkind), dimension(-num_pad+1:,-num_pad+1:,-num_pad+1:), intent(inout) :: f
    type(decomp_info), intent(in) :: gp
    if (.not. periodicBCs(1)) then
      f(0,:,:) = 2*f(1,:,:) - f(2,:,:)
      f(gp%xsz(1)+1,:,:) = 2*f(gp%xsz(1),:,:) - f(gp%xsz(1)-1,:,:)
    end if

    if (.not. periodicBCs(2)) then
      if (gp%xst(2) == 1) then
        f(:,0,:) = 2*f(:,1,:) - f(:,2,:)
      end if
      if (gp%xen(2) == gp%ysz(2)) then
        f(:,gp%xsz(2)+1,:) = 2*f(:,gp%xsz(2),:) - f(:,gp%xsz(2)-1,:)
      end if
    end if

    if (.not. periodicBCs(3)) then
      if (gp%xst(3) == 1) then
        f(:,:,0) = 2*f(:,:,1) - f(:,:,2)
      end if
      if (gp%xen(3) == gp%zsz(3)) then
        f(:,:,gp%xsz(3)+1) = 2*f(:,:,gp%xsz(3)) - f(:,:,gp%xsz(3)-1)
      end if
    end if

  end subroutine 

  subroutine advanceTime(this)
    class(enrichmentOperator), intent(inout) :: this 

  end subroutine

  subroutine renderVelocity(this)
    class(enrichmentOperator), intent(inout), target :: this
    real(rkind) :: Lx
    real(rkind), dimension(:,:), pointer :: haloBuffY, haloBuffZ 
   
    haloBuffY => null()
    haloBuffZ => null()

    Lx = xDom(2) - xDom(1)
    
    ! Zero the velocity arrays
    this%smallScales%u  = 0.d0
    this%smallScales%v  = 0.d0
    this%smallScales%wC = 0.d0
     
    ! STEP 1: Exchange Gabor modes from neighbors that have influence on your domain 
    call message(2,'Exchanging halo modes')
    call this%sendRecvHaloModes(this%modeData, 'y', this%renderModeData)
    call this%sendRecvHaloModes(this%renderModeData, 'z')

    ! Step 2: Render velocity 
      call this%renderLocalVelocity(this%renderModeData(:,1), &
        this%renderModeData(:,2),  this%renderModeData(:,3), &
        this%renderModeData(:,4),  this%renderModeData(:,5), &
        this%renderModeData(:,6),  this%renderModeData(:,7), &
        this%renderModeData(:,8),  this%renderModeData(:,9), &
        this%renderModeData(:,10), this%renderModeData(:,11), &
        this%renderModeData(:,12))

    ! Step 3: Add x-periodic contribution
    if (periodicBCs(1)) then
      call this%renderLocalVelocity(this%renderModeData(:,1) + Lx, &
        this%renderModeData(:,2),  this%renderModeData(:,3), &
        this%renderModeData(:,4),  this%renderModeData(:,5), &
        this%renderModeData(:,6),  this%renderModeData(:,7), &
        this%renderModeData(:,8),  this%renderModeData(:,9), &
        this%renderModeData(:,10), this%renderModeData(:,11), &
        this%renderModeData(:,12))
      call this%renderLocalVelocity(this%renderModeData(:,1) - Lx, &
        this%renderModeData(:,2),  this%renderModeData(:,3), &
        this%renderModeData(:,4),  this%renderModeData(:,5), &
        this%renderModeData(:,6),  this%renderModeData(:,7), &
        this%renderModeData(:,8),  this%renderModeData(:,9), &
        this%renderModeData(:,10), this%renderModeData(:,11), &
        this%renderModeData(:,12))
    end if

    ! TODO:
    ! Step 2.b: interpolate wC to w. Do we need to get uhat, vhat, what from u,
    ! v, w? 

    ! Aditya
    ! STEP 3: Impose boundary condition (no-penetration BC)
    if (this%imposeNoPenetrationBC) then
      call this%smallScales%projectToFixBC()
    end if

    ! Aditya 
    ! STEP 4: Compute pressure (in needed)
    if (this%renderPressure) then 
      call this%smallScales%computePressure()
    end if 

    nullify(haloBuffY,haloBuffZ)
  end subroutine  

  subroutine wrapupTimeStep(this)
    class(enrichmentOperator), intent(inout) :: this 

    if (mod(this%tid,this%tidRender) == 0) then 
      call tic()
      call this%renderVelocity()
      call toc('Velocity rendering took')
    end if 

    if (mod(this%tid,this%tio) == 0) then 
      call this%dumpSmallScales()
    end if 

    this%tid = this%tid + 1

    !call this%updateLargeScales(timeAdvance=.false.)

  end subroutine

  subroutine dumpSmallScales(this)
    class(enrichmentOperator), intent(inout) :: this 
    call this%smallScales%dumpFullField(this%smallScales%u,"uGab")
    call this%smallScales%dumpFullField(this%smallScales%v,"vGab")
    call this%smallScales%dumpFullField(this%smallScales%wC,"wGab")
  end subroutine

  function continueSimulation(this) result(doIcontinue)
    class(enrichmentOperator), intent(inout) :: this 
    logical :: doIContinue

    doIContinue = .false. 
    if (this%tid < this%tidStop) doIContinue = .true. 

  end function  
   
end module
