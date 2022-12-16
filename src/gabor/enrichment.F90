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
  use interpolatorMod,         only: interpolator
  implicit none

  private
  integer :: nthreads
  integer, dimension(4) :: neighbor
  integer, dimension(2) :: coords, dims
  logical, dimension(3) :: periodicBCs
  real(rkind), parameter :: tol = 1.d-13

  real(rkind), dimension(2) :: xDom, yDom, zDom

  public :: enrichmentOperator, xDom, yDom, zDom, nthreads, interpAndAddGrids 

  type :: enrichmentOperator
    !private
    real(rkind), dimension(:,:), pointer :: modeData
    real(rkind), dimension(:,:,:), allocatable :: rawModeData
    real(rkind), dimension(:), pointer, public :: kx, ky, kz
    real(rkind), dimension(:), pointer, public :: x, y, z
    real(rkind), dimension(:), pointer, public :: uhatR, uhatI, vhatR, &
      vhatI, whatR, whatI
    real(rkind), dimension(:), pointer, public :: T
    real(rkind) :: dt
    integer :: nvars
    logical :: isStratified = .false.

    ! Large scale data at mode locations
    real(rkind), dimension(:), allocatable :: KE_loc, L_loc

    ! Halo-padded arrays
    real(rkind), dimension(:,:,:), allocatable :: uh, vh, wh
    real(rkind), dimension(:,:,:), allocatable :: Th
    real(rkind), dimension(:,:,:), allocatable :: KE, L, KEh, Lh
    real(rkind), dimension(:,:,:,:), allocatable :: duidxj_h

    integer, public :: nxsupp, nysupp, nzsupp
    type(igrid), pointer :: largeScales, smallScales
    type(QHmesh), public :: QHgrid 
    real(rkind), dimension(2) :: PExbound, PEybound, PEzbound ! MPI rank boundaries
    real(rkind), dimension(:,:,:,:), pointer :: duidxj_LS
    real(rkind) :: kmin, kmax
    logical :: imposeNoPenetrationBC = .false.

    ! Initialization parameters
    integer     :: nk, ntheta, nmodes, nmodesGlobal
    real(rkind) :: scalefact, Anu, numolec, ctauGlobal
    logical     :: renderPressure = .false.
    integer     :: tidRender, tio , tidStop, tidsim
    integer, pointer :: tid

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

    integer :: activeIndex 
    contains
      procedure          :: init
      procedure          :: reinit
      procedure          :: destroy
      procedure          :: generateIsotropicModes
      procedure          :: strainModes
      procedure          :: getLargeScaleDataAtModeLocation
      procedure          :: advanceTime
      procedure          :: renderVelocity
      procedure          :: renderLocalVelocity
      procedure          :: updateLargeScaleHaloes
      procedure          :: wrapupTimeStep
      procedure          :: dumpSmallScales
      procedure          :: continueSimulation
      procedure          :: dumpData
      procedure          :: generateModes  
      procedure, private :: doDebugChecks
      procedure, private :: applyPeriodicOffsets

      ! MPI communication stuff
      procedure, private :: sendRecvHaloModes
      procedure, private :: sortAndExchangeModes 
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
    logical :: doNotRenderInitialCondition = .false.
    logical :: xPeriodic = .true., yPeriodic = .true., zPeriodic = .true.
    real(rkind) :: kminFact = 1.d0
    logical :: readGradients = .false.
    
    namelist /IO/      outputdir, writeIsotropicModes, readGradients
    namelist /GABOR/   nk, ntheta, scalefact, ctauGlobal, Anu, numolec, &
      strainInitialCondition, doNotRenderInitialCondition, &
      xPeriodic, yPeriodic, zPeriodic, kminFact
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
    this%tid => this%smallScales%step 
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
    call getDomainBoundaries(xDom,yDom,zDom,largeScales%mesh)

    this%duidxj_LS => largeScales%duidxjC
    this%isStratified = this%largeScales%isStratified
    if (this%isStratified) then
      this%nvars = 13
    else
      this%nvars = 12
    end if

    allocate(this%KE(this%smallScales%gpC%xsz(1),this%smallScales%gpC%xsz(2),&
      this%smallScales%gpC%xsz(3)))
    allocate(this%L(this%smallScales%gpC%xsz(1),this%smallScales%gpC%xsz(2),&
      this%smallScales%gpC%xsz(3)))
    call getLargeScaleParams(this%KE,this%L,this%smallScales)

    call this%QHgrid%init(inputfile, this%largeScales)

    this%nxsupp = nint(2*this%QHgrid%dx/this%smallScales%dx)
    this%nysupp = nint(2*this%QHgrid%dy/this%smallScales%dy)
    this%nzsupp = nint(2*this%QHgrid%dz/this%smallScales%dz)

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
    this%kmin = kminFact*this%kmin

    ! Allocate memory for mode data
    !allocate(this%rawData(this%nmodes,this%nvars))
    allocate(this%rawModeData(this%nmodes,this%nvars,0:1))
    this%activeIndex = 0
#include "enrichment_files/togglePointer.F90"
   
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
    this%PExbound = [this%QHgrid%xE(1), this%QHgrid%xE(this%QHgrid%gpC%xsz(1) + 1)]
    this%PEybound = [this%QHgrid%yE(1), this%QHgrid%yE(this%QHgrid%gpC%xsz(2) + 1)]
    this%PEzbound = [this%QHgrid%zE(1), this%QHgrid%zE(this%QHgrid%gpC%xsz(3) + 1)]
    call MPI_Cart_Get(DECOMP_2D_COMM_CART_X,2,dims,periodicBCs(2:3),coords,ierr)
    periodicBCs(1) = xPeriodic
    periodicBCs(2) = yPeriodic
    periodicBCs(3) = zPeriodic
    call getneighbors(neighbor,this%PEybound,this%PEzbound,periodicBCs)
   
  end subroutine

  subroutine generateModes(this) 
    class(enrichmentOperator), intent(inout) :: this 

    ! Get halo'ed large scales
    call this%updateLargeScaleHaloes(initializing=.true.) 
    
    ! Initialize the Gabor modes
    call this%generateIsotropicModes()

    if (this%writeIsotropicModes) then
      call message(1, 'Writing modes to disk.')
      call this%dumpData(this%x,this%y,this%z,this%kx,this%ky,this%kz, &
        this%uhatR,this%uhatI,this%vhatR,this%vhatI,this%whatR,this%whatI, &
        this%KE_loc, this%L_loc)
    end if
    if (this%strainInitialCondition) call this%strainModes()
  end subroutine 

  subroutine reinit(this)
    class(enrichmentOperator), intent(inout), target :: this 
    
    call getLargeScaleParams(this%KE,this%L,this%largeScales)
    !call this%updateLargeScales(timeAdvance=.false.,initializing=.true.) 

    ! Compute the number of modes
    this%nmodes = this%nk*this%ntheta * &
      this%QHgrid%gpC%xsz(1)*this%QHgrid%gpC%xsz(2)*this%QHgrid%gpC%xsz(3)
    
    ! Allocate memory for mode data
    !allocate(this%rawData(this%nmodes,this%nvars))
    if (allocated(this%rawModeData)) deallocate(this%rawModeData)
    allocate(this%rawModeData(this%nmodes,this%nvars,0:1))
    this%activeIndex = 0
#include "enrichment_files/togglePointer.F90"
   
    ! Initialize the Gabor modes
    call this%generateIsotropicModes()

    if (this%writeIsotropicModes) then
      call message(1, 'Writing modes to disk.')
      call this%dumpData(this%x,this%y,this%z,this%kx,this%ky,this%kz, &
        this%uhatR,this%uhatI,this%vhatR,this%vhatI,this%whatR,this%whatI, &
        this%KE_loc, this%L_loc)
    end if
    if (this%strainInitialCondition) call this%strainModes()
   
    !call this%wrapupTimeStep(doNotRender = doNotRenderInitialCondition)
  end subroutine

  subroutine destroy(this)
    class(enrichmentOperator), intent(inout) :: this 

    if (associated(this%largeScales)) nullify(this%largeScales)
    if (associated(this%smallScales)) nullify(this%smallScales)
    if (allocated(this%rawmodeData))    deallocate(this%rawmodeData)
    if (associated(this%modeData))    nullify(this%modeData)
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
    if (allocated(this%KEh))       deallocate(this%KEh)
    if (allocated(this%KE))       deallocate(this%KE)
    if (allocated(this%Lh))       deallocate(this%Lh)
    if (allocated(this%L))       deallocate(this%L)
    if (allocated(this%Th))       deallocate(this%Th)
    if (allocated(this%duidxj_h)) deallocate(this%duidxj_h)
    if (allocated(this%utmp))     deallocate(this%utmp)
    if (allocated(this%vtmp))     deallocate(this%vtmp)
    if (allocated(this%wtmp))     deallocate(this%wtmp)
    if (allocated(this%renderModeData)) deallocate(this%renderModeData)

    call this%QHgrid%destroy()
  end subroutine 

  subroutine updateLargeScaleHaloes(this,initializing)
    class(enrichmentOperator), intent(inout) :: this 
    integer :: i
    logical, intent(in), optional :: initializing
    logical :: init 

    init = .false. 
    if (present(initializing)) init = initializing
    ! Ryan: 
    ! STEP 1: Figure out how you want to advance the large-scales
    ! Either run a PadeOps time-step or, read in from a file 
    ! If you read in from a file, you would want to avoid doing this repeatedly. 
    ! Perhaps read in 20 snapshots at a time and so on..

    !call this%largeScales%fixGradientsForPeriodicity(periodicBCs)

    ! Aditya: 
    ! STEP 2: Generate halo'd velocities
    if(.not. allocated(this%uh))       call this%largeScales%pg%alloc_array(this%uh)
    if(.not. allocated(this%vh))       call this%largeScales%pg%alloc_array(this%vh)
    if(.not. allocated(this%wh))       call this%largeScales%pg%alloc_array(this%wh)
    if(.not. allocated(this%duidxj_h)) call this%largeScales%pg%alloc_array(this%duidxj_h,9)
    call this%largeScales%HaloUpdateVelocities(this%uh, this%vh, this%wh, &
      this%duidxj_h)

    ! The halo'ed KE and L only matter if initializing the modes
    if (init) then
      if(.not. allocated(this%KEh))      call this%smallScales%pg%alloc_array(this%KEh)
      if(.not. allocated(this%Lh))       call this%smallScales%pg%alloc_array(this%Lh)
      call this%smallScales%haloUpdateField(this%KE,this%KEh)
      call this%smallScales%haloUpdateField(this%L, this%Lh)
      call fixGhostIfNonPeriodic(this%KEh,this%smallScales%gpC,enforcePositivity=.true.)
      call fixGhostIfNonPeriodic(this%Lh,this%smallScales%gpC,enforcePositivity=.true.)
    end if

    if (this%isStratified) then
      if(.not. allocated(this%Th))       call this%largeScales%pg%alloc_array(this%Th)
      call this%largeScales%haloUpdateField(this%largeScales%T,this%Th)
      call fixGhostIfNonPeriodic(this%Th,this%largeScales%gpC)
    end if

    call fixGhostIfNonPeriodic(this%uh,this%largeScales%gpC)
    call fixGhostIfNonPeriodic(this%vh,this%largeScales%gpC)
    call fixGhostIfNonPeriodic(this%wh,this%largeScales%gpC)
    do i = 1,9
      call fixGhostIfNonPeriodic(this%duidxj_h(:,:,:,i),this%largeScales%gpC)
    end do
  end subroutine

  subroutine fixGhostIfNonPeriodic(f,gp,enforcePositivity)
    use decomp_2d, only: decomp_info
    use procgrid_mod, only: num_pad
    real(rkind), dimension(-num_pad+1:,-num_pad+1:,-num_pad+1:), intent(inout) :: f
    type(decomp_info), intent(in) :: gp
    logical, intent(in), optional :: enforcePositivity
    logical :: makePositive

    makePositive = .false.
    if (present(enforcePositivity)) makePositive = enforcePositivity

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

    if (makePositive) then
      where(f < 0.d0) f = 0.d0
    end if
  end subroutine 

  subroutine advanceTime(this)
    use GaborModeRoutines, only: rk4Step
    class(enrichmentOperator), intent(inout) :: this 
    real(rkind), dimension(3) :: k, uRtmp, uItmp, x, Ui
    real(rkind) :: L, KE
    real(rkind), dimension(3,3) :: duidxj
    integer :: n

    do n = 1,this%nmodes
      call this%getLargeScaleDataAtModeLocation(n,duidxj,Ui)
      k     = [this%kx(n)   , this%ky(n)   , this%kz(n)   ]
      uRtmp = [this%uhatR(n), this%vhatR(n), this%whatR(n)]
      uItmp = [this%uhatI(n), this%vhatI(n), this%whatI(n)]
      x     = [this%x(n)    , this%y(n)    , this%z(n)    ]
      
      call rk4Step(uRtmp,uItmp,k,x,this%dt,this%Anu,this%KE_loc(n),this%L_loc(n),this%numolec,duidxj,Ui)
      
      this%kx(n)    = k(1)
      this%ky(n)    = k(2)
      this%kz(n)    = k(3)

      this%uhatR(n) = uRtmp(1)
      this%uhatI(n) = uItmp(1)
      this%vhatR(n) = uRtmp(2)
      this%vhatI(n) = uItmp(2)
      this%whatR(n) = uRtmp(3)
      this%whatI(n) = uItmp(3)
      
      this%x(n)     = x(1)
      this%y(n)     = x(2)
      this%z(n)     = x(3)
    end do
    call this%sortAndExchangeModes(this%y, this%PEybound(1), this%PEybound(2), &
            neighbor(1), neighbor(2))    
    call this%sortAndExchangeModes(this%z, this%PEzbound(1), this%PEzbound(2), & 
            neighbor(3), neighbor(4))
   
    !if ((periodicBCs(1) .or. periodicBCs(2) .or. periodicBCs(3)) .and. &
    !  (abs(xDom(1) - this%PExbound(1)) < tol .or. abs(xDom(2) - this%PExbound(2)) < tol .or. &
    !   abs(yDom(1) - this%PEybound(1)) < tol .or. abs(yDom(2) - this%PEybound(2)) < tol .or. &
    !   abs(zDom(1) - this%PEzbound(1)) < tol .or. abs(zDom(2) - this%PEzbound(2)) < tol )) then
    !  call this%applyPeriodicOffsets()
    !end if
    if ((periodicBCs(1) .and. (abs(xDom(1) - this%PExbound(1)) < tol .or. &
                               abs(xDom(2) - this%PExbound(2)) < tol)) .or. &
        (periodicBCs(2) .and. (abs(yDom(1) - this%PEybound(1)) < tol .or. &
                               abs(yDom(2) - this%PEybound(2)) < tol)) .or. &
        (periodicBCs(3) .and. (abs(zDom(1) - this%PEzbound(1)) < tol .or. &
                               abs(zDom(2) - this%PEzbound(2)) < tol)) ) then
      call this%applyPeriodicOffsets()
    end if

    if (this%debugChecks) call this%doDebugChecks()
  end subroutine

  subroutine applyPeriodicOffsets(this)
    class(enrichmentOperator), intent(inout) :: this
    real(rkind) :: Lx, Ly, Lz
    integer :: n

    Lx = xDom(2) - xDom(1)
    Ly = yDom(2) - yDom(1)
    Lz = zDom(2) - zDom(1)

    do n = 1,this%nmodes
      ! X
      if (this%x(n) < xDom(1)) then
        call assert(periodicBCs(1),'periodicBCs(1)')
        call assert(abs(this%QHgrid%xE(this%QHgrid%gpC%xsz(1)+1) - xDom(2)) < tol,&
          'abs(this%QHgrid%xE(this%QHgrid%gpC%xsz(1)+1) - xDom(2)) < tol')
        this%x(n) = this%x(n) + Lx
      end if

      if (this%x(n) > xDom(2)) then
        call assert(periodicBCs(1),'periodicBCs(1)')
        call assert(abs(this%QHgrid%xE(1) - xDom(1)) < tol,&
          'abs(this%QHgrid%xE(1) - xDom(1)) < tol')
        this%x(n) = this%x(n) - Lx
      end if
      
      ! Y
      if (this%y(n) < yDom(1)) then
        call assert(periodicBCs(2),'periodicBCs(2)')
        call assert(abs(this%QHgrid%yE(this%QHgrid%gpC%xsz(2)+1) - yDom(2)) < tol,&
          'abs(this%QHgrid%yE(this%QHgrid%gpC%xsz(2)+1) - yDom(2)) < tol')
        this%y(n) = this%y(n) + Ly
      end if

      if (this%y(n) > yDom(2)) then
        call assert(periodicBCs(2),'periodicBCs(2)')
        call assert(abs(this%QHgrid%yE(1) - yDom(1)) < tol,&
          'abs(this%QHgrid%yE(1) - yDom(1)) < tol')
        this%y(n) = this%y(n) - Ly
      end if
      
      ! Z
      if (this%z(n) < zDom(1)) then
        call assert(periodicBCs(3),'periodicBCs(3)')
        call assert(abs(this%QHgrid%zE(this%QHgrid%gpC%xsz(3)+1) - zDom(2)) < tol,&
          'abs(this%QHgrid%zE(this%QHgrid%gpC%xsz(3)+1) - zDom(2)) < tol')
        this%z(n) = this%z(n) + Lz
      end if

      if (this%z(n) > zDom(2)) then
        call assert(periodicBCs(3),'periodicBCs(3)')
        call assert(abs(this%QHgrid%zE(1) - zDom(1)) < tol,&
          'abs(this%QHgrid%zE(1) - zDom(1)) < tol')
        this%z(n) = this%z(n) - Lz
      end if
      
    end do
  end subroutine

  subroutine sortAndExchangeModes(this, coor, coorMin, coorMax, neighLo, neighHi)
    class(enrichmentOperator), intent(inout), target :: this
    integer :: n, inactiveIndex, iterSelf, iterLo, iterHi
    integer :: howmanyLo, howmanyHi
    real(rkind), dimension(:), intent(in) :: coor
    real(rkind), intent(in) :: coorMin, coorMax
    integer, intent(in) :: neighLo, neighHi
    real(rkind), dimension(:,:), allocatable :: sendBufferLo, sendBufferHi
    real(rkind), dimension(:,:), allocatable :: recvBufferLo, recvBufferHi
    real(rkind), dimension(:,:,:), allocatable :: tmp
    integer, dimension(4) :: sendReq, recvReq, requests
    integer :: tag, ierr, sz1

    inactiveIndex = mod(this%activeIndex+1,2)
    
    allocate(sendBufferLo(this%nmodes,this%nvars))
    allocate(sendBufferHi(this%nmodes,this%nvars))

    iterSelf = 0
    iterLo = 0
    iterHi = 0

    do n = 1,this%nModes
        
        if ((coor(n) <= coorMax) .and. (coor(n) > coorMin)) then 
            iterSelf = iterSelf + 1
            this%rawModedata(iterSelf,:,inactiveIndex) = this%rawModedata(n,:,this%activeIndex) 
        else
            if (coor(n) > coorMax) then 
                iterHi = iterHi + 1
                sendBufferHi(iterHi,:) = this%rawModedata(n,:,this%activeIndex) 

            else if (coor(n) <= coorMin) then 
                iterLo = iterLo + 1
                sendBufferLo(iterLo,:) = this%rawModedata(n,:,this%activeIndex) 
            
            else
                print*, coor(n), coorMax, coorMin
                print*, "This should not happen"
                stop 
            end if 
        end if 

    end do


    ! isend howmany to hi
    ! isend howmany to lo
    tag = 0 
    ! Send to lower rank
    call MPI_ISend(iterLo,1,MPI_INTEGER,neighLo,tag,&
      DECOMP_2D_COMM_CART_X,sendReq(1),ierr)
    ! Send to higher rank
    call MPI_ISend(iterHi,1,MPI_INTEGER,neighHi,tag,&
      DECOMP_2D_COMM_CART_X,sendReq(2),ierr)
    
    ! irecv howmany from lo
    ! irecv howmany from hi
    call MPI_IRecv(howManyLo,1,MPI_INTEGER,neighLo,tag,&
      DECOMP_2D_COMM_CART_X,recvReq(1),ierr)
    ! Receive from higher rank
    call MPI_IRecv(howManyHi,1,MPI_INTEGER,neighHi,tag,&
      DECOMP_2D_COMM_CART_X,recvReq(2),ierr)
    
    ! isend data to lo
    ! isend data to hi 
    call MPI_ISend(sendBufferLo(1:iterLo,:),this%nvars*iterLo,mpirkind,neighLo,tag,&
        DECOMP_2D_COMM_CART_X,sendReq(3),ierr) 
    call MPI_ISend(sendBufferHi(1:iterHi,:),this%nvars*iterHi,mpirkind,neighHi,tag,&
        DECOMP_2D_COMM_CART_X,sendReq(4),ierr) 
    
    requests(1:2) = recvReq(1:2)
    requests(3:4) = sendReq(1:2)
    call MPI_WaitAll(4,requests,MPI_STATUSES_IGNORE,ierr)

    allocate(recvBufferLo(howmanyLo,this%nvars))
    allocate(recvBufferHi(howmanyHi,this%nvars))

    ! irecv data from lo
    ! irecv data from hi 
    call MPI_IRecv(recvBufferLo,size(recvBufferLo),mpirkind,neighLo,tag,&
      DECOMP_2D_COMM_CART_X,recvReq(3),ierr)
    call MPI_IRecv(recvBufferHi,size(recvBufferHi),mpirkind,neighHi,tag,&
      DECOMP_2D_COMM_CART_X,recvReq(4),ierr)
  
    requests(1:2) = recvReq(3:4)
    requests(3:4) = sendReq(3:4)
    call MPI_WaitAll(4,requests, MPI_STATUSES_IGNORE,ierr)

    if (size(this%rawModedata,1) < iterSelf+howmanyLo+howmanyHi) then
        sz1 = size(this%rawModeData,1)
        
        ! Temporarily copy rawModeData to tmp array
        allocate(tmp(sz1,size(this%rawModeData,2),size(this%rawModeData,3)))
        tmp = this%rawModeData

        ! Deallocate and reallocate rawModeData
        deallocate(this%rawModeData)
        allocate(this%rawModeData(iterSelf+howManyLo+howManyHi,this%nvars,0:1))

        ! Copy tmp data back to rawModeData
        this%rawModeData(1:sz1,:,:) = tmp

        ! Clear memory
        deallocate(tmp)
    end if 

    this%rawModedata(iterSelf+1:iterSelf+howmanyLo,:,inactiveIndex) = recvBufferLo 
    this%rawModedata(iterSelf+howmanyLo+1:iterSelf+howmanyLo+howmanyHi,:,inactiveIndex) = recvBufferHi 

    this%nModes = iterSelf+howmanyLo+howmanyHi

#include "enrichment_files/togglePointer.F90"

    deallocate(sendBufferLo, sendBufferHi, recvBufferLo, recvBufferHi)

  end subroutine 

  subroutine renderVelocity(this)
    class(enrichmentOperator), intent(inout), target :: this
    real(rkind) :: Lx
    real(rkind), dimension(:,:), pointer :: haloBuffY, haloBuffZ 
    
    call message(1,"Rendering the Gabor-induced velocity field")
   
    haloBuffY => null()
    haloBuffZ => null()

    Lx = xDom(2) - xDom(1)
    
    ! Zero the velocity arrays
    this%smallScales%u  = 0.d0
    this%smallScales%v  = 0.d0
    this%smallScales%wC = 0.d0
     
    ! STEP 1: Exchange Gabor modes from neighbors that have influence on your domain 
    if (allocated(this%renderModeData)) deallocate(this%renderModedata)
    allocate(this%renderModeData(this%nmodes, size(this%ModeData,2)))
    !allocate(this%renderModeData(size(this%modeData,1), size(this%ModeData,2)))
    this%renderModeData = this%modeData
    call this%sendRecvHaloModes(this%renderModeData, 'y')
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

  subroutine wrapupTimeStep(this,doNotRender)
    class(enrichmentOperator), intent(inout) :: this
    logical :: noRender
    logical, intent(in), optional :: doNotRender

    noRender = .false.
    if (present(doNotRender)) noRender = doNotRender
    if (noRender) then
      continue
    else
      if (mod(this%tid,this%tidRender) == 0 .or. this%tid == this%tidStop-1) then
          call tic()
          call this%renderVelocity()
          call toc('Velocity rendering took')
      end if 

      if (mod(this%tid,this%tio) == 0 .or. this%tid == this%tidStop-1) then
        call this%dumpSmallScales()
      end if
    end if 
    
    this%tid = this%tid + 1


  end subroutine

  subroutine dumpSmallScales(this)
    class(enrichmentOperator), intent(inout) :: this 
    call this%smallScales%dumpFullField(this%smallScales%u,"uVel")
    call this%smallScales%dumpFullField(this%smallScales%v,"vVel")
    call this%smallScales%dumpFullField(this%smallScales%wC,"wVel")
  end subroutine

  function continueSimulation(this) result(doIcontinue)
    class(enrichmentOperator), intent(inout) :: this 
    logical :: doIContinue

    doIContinue = .false. 
    if (this%tid < this%tidStop) doIContinue = .true. 

  end function  

  subroutine interpAndAddGrids(sourceGrid, destGrid, interp)
    class(igrid), intent(in) :: sourceGrid
    class(igrid), intent(inout) :: destGrid
    class(interpolator) :: interp 
    integer :: iter 
     
    call interp%LinInterp3D(sourceGrid%u ,destGrid%rbuffxC(:,:,:,1))
    destGrid%u = destGrid%u + destGrid%rbuffxC(:,:,:,1)

    call interp%LinInterp3D(sourceGrid%v ,destGrid%rbuffxC(:,:,:,1))
    destGrid%v = destGrid%v + destGrid%rbuffxC(:,:,:,1)

    call interp%LinInterp3D(sourceGrid%wC ,destGrid%rbuffxC(:,:,:,1))
    destGrid%wC = destGrid%wC + destGrid%rbuffxC(:,:,:,1)

    do iter = 1,9
        call interp%LinInterp3D(sourceGrid%duidxjC(:,:,:,iter) ,destGrid%rbuffxC(:,:,:,1))
        destGrid%duidxjC(:,:,:,iter) = destGrid%duidxjC(:,:,:,iter) + destGrid%rbuffxC(:,:,:,1)
    end do 

  end subroutine 

end module
