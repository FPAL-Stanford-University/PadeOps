module enrichmentMod
  use kind_parameters,    only: rkind, clen
  use incompressibleGrid, only: igrid
  use QHmeshMod,          only: QHmesh
  use exits,              only: message
  use constants,          only: pi
  use fortran_assert,     only: assert
  use omp_lib
  use gaborHooks
  implicit none

  integer :: nthreads

  type :: enrichmentOperator
    real(rkind), dimension(:), allocatable :: kx, ky, kz
    real(rkind), dimension(:), allocatable :: x, y, z
    real(rkind), dimension(:), allocatable :: uhatR, uhatI, vhatR, vhatI, whatR, whatI
    real(rkind) :: dt

    ! Halo-padded arrays
    real(rkind), dimension(:,:,:), allocatable :: uh, vh, wh
    real(rkind), dimension(:,:,:,:), allocatable :: duidxj_h

    integer :: nxsupp, nysupp, nzsupp
    type(igrid), pointer :: largeScales, smallScales
    type(QHmesh) :: QHgrid 
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
    real(rkind), dimension(:,:,:,:), allocatable :: utmp,vtmp,wtmp

    ! Misc
    logical :: debugChecks = .false.

    contains
      procedure          :: init
      procedure          :: destroy
      procedure          :: generateIsotropicModes
      procedure          :: strainModes
      procedure          :: getLargeScaleDataAtModeLocation
      procedure          :: advanceTime
      procedure          :: renderVelocity
      procedure, private :: renderLocalVelocity
      procedure          :: updateLargeScales
      procedure          :: wrapupTimeStep
      procedure          :: continueSimulation
      procedure          :: dumpData
      procedure, private :: doDebugChecks
  end type

contains
  include 'enrichment_files/initializeGaborModes.F90'
  include 'enrichment_files/renderVelocity.F90'
  include 'enrichment_files/gaborIO.F90'
  include 'enrichment_files/debugChecks.F90'

  subroutine init(this,smallScales,largeScales,inputfile,Lx,Ly,Lz)
    use GaborModeRoutines, only: computeKminKmax
    
    class(enrichmentOperator), intent(inout) :: this 
    class(igrid), intent(inout), target :: smallScales, largeScales
    real(rkind), intent(in) :: Lx, Ly, Lz
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
    
    namelist /IO/      outputdir, writeIsotropicModes
    namelist /GABOR/   nk, ntheta, scalefact, ctauGlobal, Anu, numolec
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

    ! IO
    this%outputdir = outputdir
    this%writeIsotropicModes = writeIsotropicModes

    ! STEP  0: Get filenames for u, v and w from (input file?) for initialization
    ! Use the same file-naming convention for PadeOps restart files: 
    ! RESTART_Run00_u.000000
    call this%largeScales%initLargeScales(tidInit,this%largeScales%runID)
    this%duidxj_LS => largeScales%duidxjC

    ! Get halo-padded velocity arrays
    call this%updateLargeScales(timeAdvance=.false.) 

    ! Ryan 
    ! STEP 1: Finish the QH region code to fill Gabor Modes (kx, ky, kz, x, y, z, uhat, ...)
    call this%QHgrid%init(inputfile,this%largeScales)
    call getLargeScaleParams(this%QHgrid%KE,this%QHgrid%L,this%largeScales)

    this%nxsupp = 2 * this%smallScales%nx/this%largeScales%nx * &
      nint(this%QHgrid%dx / this%largeScales%dx)
    this%nysupp = 2 * this%smallScales%ny/this%largeScales%ny * &
      nint(this%QHgrid%dy / this%largeScales%dy)
    this%nzsupp = 2 * this%smallScales%nz/this%largeScales%nz * &
      nint(this%QHgrid%dz / this%largeScales%dz)

    this%nmodes = this%nk*this%ntheta * &
      this%QHgrid%gpC%xsz(1)*this%QHgrid%gpC%xsz(2)*this%QHgrid%gpC%xsz(3)
    this%nmodesGlobal = this%nk*this%ntheta * &
      this%QHgrid%nx*this%QHgrid%ny*this%QHgrid%nz
    
    ! Compute kmin and kmax based on LES and high-resolution grids 
    call computeKminKmax(Lx, Ly, Lz, &
      this%largeScales%nx, this%largeScales%ny, this%largeScales%nz, &
      this%smallScales%nx, this%smallScales%ny, this%smallScales%nz, &
      this%kmin, this%kmax)

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
    
    ! Initialize the Gabor modes
    call this%generateIsotropicModes()
    if (this%writeIsotropicModes) call this%dumpData()
    call this%strainModes()

    call this%wrapupTimeStep()

  end subroutine

  subroutine destroy(this)
    class(enrichmentOperator), intent(inout) :: this 

    if (associated(this%largeScales)) nullify(this%largeScales)
    if (associated(this%smallScales)) nullify(this%smallScales)
    if (allocated(this%uhatR))    deallocate(this%uhatR)
    if (allocated(this%uhatI))    deallocate(this%uhatI)
    if (allocated(this%vhatR))    deallocate(this%vhatR)
    if (allocated(this%vhatI))    deallocate(this%vhatI)
    if (allocated(this%whatR))    deallocate(this%whatR)
    if (allocated(this%whatI))    deallocate(this%whatI)
    if (allocated(this%kx))       deallocate(this%kx)
    if (allocated(this%ky))       deallocate(this%ky)
    if (allocated(this%kz))       deallocate(this%kz)
    if (allocated(this%x))        deallocate(this%x)
    if (allocated(this%y))        deallocate(this%y)
    if (allocated(this%z))        deallocate(this%z)
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

    ! Ryan: 
    ! STEP 1: Figure out how you want to advance the large-scales
    ! Either run a PadeOps time-step or, read in from a file 
    ! If you read in from a file, you would want to avoid doing this repeatedly. 
    ! Perhaps read in 20 snapshots at a time and so on..
    if (timeAdvance) call this%largeScales%timeAdvance(this%dt)

    ! Aditya: 
    ! STEP 2: Generate halo'd velocities
    call this%largeScales%HaloUpdateVelocities(this%uh, this%vh, this%wh, &
      this%duidxj_h)

  end subroutine 

  subroutine advanceTime(this)
    class(enrichmentOperator), intent(inout) :: this 

  end subroutine

  subroutine renderVelocity(this)
    class(enrichmentOperator), intent(inout) :: this 

    ! Aditya & Ryan 
    ! STEP 1: Exchange Gabor modes from neighbors that have influence on your domain 


    ! Ryan
    ! STEP 2: Render local velocity
    call this%renderLocalVelocity()

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

  end subroutine  

  subroutine wrapupTimeStep(this)
    class(enrichmentOperator), intent(inout) :: this 

    if (mod(this%tid,this%tidRender) == 0) then 
      call this%renderVelocity()
    end if 

    if (mod(this%tid,this%tio) == 0) then 
      call this%smallScales%dumpFullField(this%smallScales%u,"uGab")
      call this%smallScales%dumpFullField(this%smallScales%v,"vGab")
      call this%smallScales%dumpFullField(this%smallScales%wC,"wGab")
    end if 

    this%tid = this%tid + 1

    !call this%updateLargeScales(timeAdvance=.false.)

  end subroutine

  function continueSimulation(this) result(doIcontinue)
    class(enrichmentOperator), intent(inout) :: this 
    logical :: doIContinue

    doIContinue = .false. 
    if (this%tid < this%tidStop) doIContinue = .true. 

  end function  
      
end module
