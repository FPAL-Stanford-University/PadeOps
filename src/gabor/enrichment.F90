module enrichmentMod
  use kind_parameters,    only: rkind
  use incompressibleGrid, only: igrid
  implicit none

  type :: enrichmentOperator
    real(rkind), dimension(:),     allocatable :: kx, ky, kz
    real(rkind), dimension(:),     allocatable :: x, y, z
    real(rkind), dimension(:),     allocatable :: uhatR, uhatI, vhatR, vhatI, whatR, whatI
    real(rkind), dimension(:,:,:), allocatable :: uh, vh, wh
    type(igrid), pointer :: largeScales, smallScales 

    contains
      procedure :: init
      procedure :: destroy
      !procedure :: generateIsotropicModes
      !procedure :: strainModes
      !procedure :: timeAdvance
      !procedure :: renderVelocity
      !procedure :: updateLargeScales
    end type

contains
  subroutine init(this,smallScales,largeScales,inputfile)
    class(enrichmentOperator), intent(inout) :: this 
    class(igrid), intent(inout), target :: smallScales, largeScales
    character(len=*), intent(in) :: inputfile   
    integer :: ierr, ioUnit
    integer :: nk, ntheta
    real(rkind) :: scalefact, Anu, numolec, ctau
     
    namelist /GABOR/ nk, ntheta, scalefact, ctau, Anu, numolec

    ! Read inputfile
    ioUnit = 1
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED', iostat=ierr)
    read(unit=ioUnit, NML=GABOR)
    close(ioUnit)

    this%largeScales => largeScales 
    this%smallScales => smallScales

  end subroutine

  subroutine destroy(this)
    class(enrichmentOperator), intent(inout) :: this 

    if (associated(this%largeScales)) nullify(this%largeScales)
    if (associated(this%smallScales)) nullify(this%smallScales)
  end subroutine  
end module
