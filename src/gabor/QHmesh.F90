module QHmeshMod
  use kind_parameters,    only: rkind, clen
  use decomp_2d,          only: nproc, nrank, decomp_info, decomp_info_init, &
                                decomp_info_finalize
  use incompressibleGrid, only: igrid
  use exits,              only: message
  use reductions,         only: p_minval, p_sum
  use fortran_assert,     only: assert
  implicit none

  type QHmesh
    real(rkind), dimension(:), allocatable :: xE, yE, zE
    real(rkind), dimension(:), allocatable :: xC, yC, zC
    real(rkind) :: dx, dy, dz 
    integer :: nx, ny, nz
    integer :: isz, jsz, ksz
    integer :: ist, ien, jst, jen, kst, ken
    type(igrid), pointer :: LES
    integer, dimension(:,:,:), allocatable :: gID
    
    contains
      procedure :: init
      procedure :: destroy
  end type

  contains

    subroutine init(this,inputfile,largeScales,smallScales,xDom,yDom,zDom)
      use QHmeshRoutines, only: getMeshEdgeValues
      class(QHmesh), intent(inout) :: this
      character(len=clen), intent(in) :: inputfile
      class(igrid), intent(in)  :: smallScales, largeScales
      real(rkind), dimension(2), intent(in) :: xDom, yDom, zDom
      class(decomp_info), allocatable :: gpC
      integer :: nxLESperQH, nyLESperQH, nzLESperQH, nx, ny, nz
      integer :: ioUnit, ierr, ntotal
      integer :: i, j, k, ig, jg, kg
      real(rkind) :: dx, dy, dz

      namelist /QHMESH/ nxLESperQH, nyLESperQH, nzLESperQH

      ioUnit = 1
      open(unit=ioUnit, file=trim(inputfile), form='FORMATTED', iostat=ierr)
      read(unit=ioUnit, NML=QHMESH)
      close(ioUnit)

      call assert(mod(largeScales%nx,nxLESperQH) == 0,'nx in LES must be a '&
        //'multiple of nxLESperQH -- QHmesh.F90')      
      call assert(mod(largeScales%ny,nyLESperQH) == 0,'ny in LES must be a '&
        //'multiple of nyLESperQH -- QHmesh.F90')      
      call assert(mod(largeScales%nz,nzLESperQH) == 0,'nz in LES must be a '&
        //'multiple of nzLESperQH -- QHmesh.F90')      
      
      this%nx = largeScales%nx/nxLESperQH
      this%ny = largeScales%ny/nyLESperQH
      this%nz = largeScales%nz/nzLESperQH
      ntotal = this%nx*this%ny*this%nz

      this%dx = real(nxLESperQH,rkind)*largeScales%dx
      this%dy = real(nyLESperQH,rkind)*largeScales%dy
      this%dz = real(nzLESperQH,rkind)*largeScales%dz

      call message('*** NOTE: QHmesh class assumes x-pencil decomposition ***')

      allocate(gpC)
      call decomp_info_init(this%nx, this%ny, this%nz, gpC)

      ! Get the local domain bounds based on the small scales mesh
      dx = smallScales%mesh(2,1,1,1) - smallScales%mesh(1,1,1,1)
      dy = smallScales%mesh(1,2,1,2) - smallScales%mesh(1,1,1,2)
      dz = smallScales%mesh(1,1,2,3) - smallScales%mesh(1,1,1,3)

      this%isz = gpC%xsz(1)
      this%jsz = gpC%xsz(2)
      this%ksz = gpC%xsz(3)

      this%ist = gpC%xst(1)
      this%jst = gpC%xst(2)
      this%kst = gpC%xst(3)
      
      this%ien = gpC%xen(1)
      this%jen = gpC%xen(2)
      this%ken = gpC%xen(3)

      print*, "Number of QH regions:", this%nx, this%ny, this%nz

      ! Allocate memory
      allocate(this%xC(this%isz), this%yC(this%jsz), this%zC(this%ksz))
      allocate(this%xE(this%isz+1), this%yE(this%jsz+1), this%zE(this%ksz+1))
      allocate(this%gID(this%isz,this%jsz,this%ksz))

      ! Define the mesh        
      this%xE = getMeshEdgeValues(this%ist, this%isz+1, this%dx, xDom(1))
      this%yE = getMeshEdgeValues(this%jst, this%jsz+1, this%dy, yDom(1))
      this%zE = getMeshEdgeValues(this%kst, this%ksz+1, this%dz, zDom(1))
      
      this%xC = 0.5d0*(this%xE(2:this%isz+1)+this%xE(1:this%isz))
      this%yC = 0.5d0*(this%yE(2:this%jsz+1)+this%yE(1:this%jsz))
      this%zC = 0.5d0*(this%zE(2:this%ksz+1)+this%zE(1:this%ksz))

      this%dx = this%xE(2) - this%xE(1)
      this%dy = this%yE(2) - this%yE(1)
      this%dz = this%zE(2) - this%zE(1)

      do k = 1,this%ksz
        kg = this%kst + k - 1
        do j = 1,this%jsz
          jg = this%jst + j - 1
          do i = 1,this%isz
            ig = i
            this%gID(i,j,k) = (kg - 1)*this%ny*this%nx + (jg - 1)*this%nx + ig
          end do
        end do
      end do

      print*, nrank, "sizes:", this%isz, this%jsz, this%ksz
      print*, nrank, "starts:", this%ist, this%jst, this%kst
      print*, nrank, "max(gID):", maxval(this%gID)
      call decomp_info_finalize(gpC)
      deallocate(gpC)
    end subroutine
    
    subroutine destroy(this)
      class(QHmesh), intent(inout) :: this

      if (allocated(this%xE)) deallocate(this%xE)
      if (allocated(this%yE)) deallocate(this%yE)
      if (allocated(this%zE)) deallocate(this%zE)
      if (allocated(this%xC)) deallocate(this%xC)
      if (allocated(this%yC)) deallocate(this%yC)
      if (allocated(this%zC)) deallocate(this%zC)
      if (allocated(this%gID)) deallocate(this%gID)
    end subroutine
    
end module
