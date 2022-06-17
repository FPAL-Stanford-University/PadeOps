module QHmeshMod
  use kind_parameters,    only: rkind, clen
  use decomp_2d
  use incompressibleGrid, only: igrid
  use exits,              only: message
  implicit none

  type QHmesh
    type(decomp_info) :: gpC
    real(rkind), dimension(:), allocatable :: xE, yE, zE
    real(rkind), dimension(:), allocatable :: xC, yC, zC
    real(rkind) :: dx, dy, dz 
    integer :: nx, ny, nz
    type(igrid), pointer :: LES
    real(rkind), dimension(:,:,:), allocatable :: KE, L 
    
    contains
      procedure :: init
      procedure :: destroy
  end type

  contains

    subroutine init(this,inputfile,largeScales)
      use QHmeshRoutines, only: getMeshEdgeValues
      class(QHmesh), intent(inout) :: this
      character(len=clen), intent(in) :: inputfile
      class(igrid), intent(in), target :: largeScales
      integer :: ist, ien, jst, jen, kst, ken, isz, jsz, ksz
      integer :: nxLESperQH, nyLESperQH, nzLESperQH
      integer :: ioUnit, ierr

      namelist /QHMESH/ nxLESperQH, nyLESperQH, nzLESperQH

      ioUnit = 1
      open(unit=ioUnit, file=trim(inputfile), form='FORMATTED', iostat=ierr)
      read(unit=ioUnit, NML=QHMESH)
      close(ioUnit)

      this%LES => largeScales
      
      this%nx = this%LES%nx/nxLESperQH
      this%ny = this%LES%ny/nyLESperQH
      this%nz = this%LES%nz/nzLESperQH

      this%dx = real(nxLESperQH,rkind)*this%LES%dx
      this%dy = real(nyLESperQH,rkind)*this%LES%dy
      this%dz = real(nzLESperQH,rkind)*this%LES%dz

      call decomp_info_init(this%nx,this%ny,this%nz,this%gpC)
      
      call message('WARNING: QHmesh class assumes x-pencil decomposition')
      ist = this%gpC%xst(1)
      ien = this%gpC%xen(1)
      jst = this%gpC%xst(2)
      jen = this%gpC%xen(2)
      kst = this%gpC%xst(3)
      ken = this%gpC%xen(3)

      isz = this%gpC%xsz(1)
      jsz = this%gpC%xsz(2)
      ksz = this%gpC%xsz(3)

      ! Allocate memory
      allocate(this%xC(isz), this%yC(jsz), this%zC(ksz))
      allocate(this%xE(isz+1), this%yE(jsz+1), this%zE(ksz+1))

      allocate(this%KE(isz,jsz,ksz), this%L(isz,jsz,ksz))

      ! Define the mesh           
      this%xE = getMeshEdgeValues(ist,ien+1,this%dx)
      this%yE = getMeshEdgeValues(jst,jen+1,this%dy)
      this%zE = getMeshEdgeValues(kst,ken+1,this%dz)
      
      this%xC = 0.5d0*(this%xE(2:isz+1)+this%xE(1:isz))
      this%yC = 0.5d0*(this%yE(2:jsz+1)+this%yE(1:jsz))
      this%zC = 0.5d0*(this%zE(2:ksz+1)+this%zE(1:ksz))

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
      if (allocated(this%KE)) deallocate(this%KE)
      if (allocated(this%L))  deallocate(this%L)
      call decomp_info_finalize(this%gpC)
    end subroutine
    
end module
