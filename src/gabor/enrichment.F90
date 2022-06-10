module enrichmentMod
  use kind_parameters,    only: rkind
  use incompressibleGrid, only: igrid
  use QHmeshMod,          only: QHmesh
  implicit none

  type :: enrichmentOperator
    real(rkind), dimension(:),     allocatable :: kx, ky, kz
    real(rkind), dimension(:),     allocatable :: x, y, z
    real(rkind), dimension(:),     allocatable :: uhatR, uhatI, vhatR, vhatI, whatR, whatI
    real(rkind), dimension(:,:,:), allocatable :: uh, vh, wh
    type(igrid), pointer :: largeScales, smallScales
    integer :: tidRender, tio , tid, tidStop 

    contains
      procedure :: init
      procedure :: destroy
      procedure :: generateIsotropicModes
      !procedure :: strainModes
      procedure :: advanceTime
      procedure :: renderVelocity
      procedure :: updateLargeScales
      procedure :: wrapupTimeStep
      procedure :: continueSimulation 
  end type

contains
  subroutine init(this,smallScales,largeScales,inputfile)
    class(enrichmentOperator), intent(inout) :: this 
    class(igrid), intent(inout), target :: smallScales, largeScales
    character(len=*), intent(in) :: inputfile   
    integer :: ierr, ioUnit
    integer :: nk, ntheta
    integer :: tidRender, tio, tidStop 
    real(rkind) :: scalefact, Anu, numolec, ctau
     
    namelist /GABOR/ nk, ntheta, scalefact, ctau, Anu, numolec
    namelist /CONTROL/ tidRender, tio, tidStop 

    ! Read inputfile
    ioUnit = 1
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED', iostat=ierr)
    read(unit=ioUnit, NML=GABOR)
    read(unit=ioUnit, NML=CONTROL)
    close(ioUnit)

    this%largeScales => largeScales 
    this%smallScales => smallScales

    this%tidRender = tidRender 
    this%tio = tio
    this%tidStop = tidStop 
  end subroutine

  subroutine destroy(this)
    class(enrichmentOperator), intent(inout) :: this 

    if (associated(this%largeScales)) nullify(this%largeScales)
    if (associated(this%smallScales)) nullify(this%smallScales)
  end subroutine 

  subroutine generateIsotropicModes(this,QHregions)
    class(enrichmentOperator), intent(inout) :: this
    type(QHmesh), intent(in) :: QHregions
    integer :: i, j, k
    real(rkind) :: xmin, ymin, zmin

    do k = 1,QHregions%nz
      zmin = QHregions%xC(k)
      do j = 1,QHregions%ny
        ymin = QHregions%yC(j)
        do i = 1,QHregions%nx
          xmin = QHregions%xC(i)

          ! Initialize modes

        end do
      end do
    end do
  end subroutine 

  subroutine updateLargeScales(this)
    class(enrichmentOperator), intent(inout) :: this 


  end subroutine 

  subroutine advanceTime(this)
    class(enrichmentOperator), intent(inout) :: this 

  end subroutine

  subroutine renderVelocity(this)
    class(enrichmentOperator), intent(inout) :: this 


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

  end subroutine

  function continueSimulation(this) result(doIcontinue)
    class(enrichmentOperator), intent(inout) :: this 
    logical :: doIContinue

    doIContinue = .false. 
    if (this%tid < this%tidStop) doIContinue = .true. 

  end function  
end module
