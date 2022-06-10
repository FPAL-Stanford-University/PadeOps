module QHmeshMod
  use kind_parameters,    only: rkind, clen
  use decomp_2d
  use incompressibleGrid, only: igrid
  implicit none

  type QHmesh
    type(decomp_info) :: gpC
    real(rkind), dimension(:), allocatable :: xE, yE, zE
    real(rkind), dimension(:), allocatable :: xC, yC, zC
    real(rkind) :: dx, dy, dz 
    integer :: nx, ny, nz
    type(igrid), pointer :: LES
    
    contains
      procedure :: init
      procedure :: destroy
  end type

  contains

    subroutine init(this,inputfile,largeScales)
      class(QHmesh), intent(inout) :: this
      character(len=clen), intent(in) :: inputfile
      class(igrid), intent(in), target :: largeScales
      integer :: nxLESperQH, nyLESperQH, nzLESperQH
      integer :: ioUnit, ierr

      namelist /QHMESH/ nxLESperQH, nyLESperQH, nzLESperQH

      open(unit=ioUnit, file=trim(inputfile), form='FORMATTED', iostat=ierr)
      read(unit=ioUnit, NML=QHMESH)
      close(ioUnit)

      this%LES => largeScales
      
      this%nx = this%LES%nx/nxLESperQH
      this%ny = this%LES%ny/nyLESperQH
      this%nz = this%LES%nz/nzLESperQH

      call decomp_info_init(this%nx,this%ny,this%nz,this%gpC)

      allocate(this%xC(this%gpC%xsz(1)),   this%yC(this%gpC%xsz(2)),   &
        this%zC(this%gpC%xsz(3)))
      allocate(this%xE(this%gpC%xsz(1)+1), this%yE(this%gpC%xsz(2)+1), &
        this%zE(this%gpC%xsz(3)+1))

    end subroutine
    
    subroutine destroy(this)
      class(QHmesh), intent(inout) :: this

      if (associated(this%LES)) nullify(this%LES)
      if (allocated(this%xE)) deallocate(this%xE)
      if (allocated(this%yE)) deallocate(this%yE)
      if (allocated(this%zE)) deallocate(this%zE)
      if (allocated(this%xC)) deallocate(this%xC)
      if (allocated(this%yC)) deallocate(this%yC)
      if (allocated(this%zC)) deallocate(this%zC)
      call decomp_info_finalize(this%gpC)
    end subroutine
end module
