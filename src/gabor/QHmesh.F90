module QHmeshMod
  use kind_parameters,    only: rkind, clen
  use decomp_2d!,          only: nproc, nrank
  use incompressibleGrid, only: igrid
  use exits,              only: message
  use reductions,         only: p_minval, p_sum
  use fortran_assert,     only: assert
  use mpi
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
      use gridtools, only: linspace
      class(QHmesh), intent(inout) :: this
      character(len=clen), intent(in) :: inputfile
      class(igrid), intent(in)  :: smallScales, largeScales
      real(rkind), dimension(2), intent(in) :: xDom, yDom, zDom
      integer :: nxLESperQH, nyLESperQH, nzLESperQH, nx, ny, nz
      integer :: ioUnit, ierr, ntotal
      integer :: i, j, k, ig, jg, kg
      integer, dimension(nproc) :: iszAll, jszAll, kszAll
      real(rkind) :: dx, dy, dz
      real(rkind) :: xst, xen, yst, yen, zst, zen
      class(decomp_info), allocatable :: gpC

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

!      allocate(gpC)
!      call decomp_info_init(this%nx,this%ny,this%nz,gpC)
      
      call message('*** NOTE: QHmesh class assumes x-pencil decomposition ***')
!      this%ist = gpC%xst(1)
!      this%ien = gpC%xen(1)
!      this%jst = gpC%xst(2)
!      this%jen = gpC%xen(2)
!      this%kst = gpC%xst(3)
!      this%ken = gpC%xen(3)
!
!      this%isz = gpC%xsz(1)
!      this%jsz = gpC%xsz(2)
!      this%ksz = gpC%xsz(3)

      ! Get the local domain bounds based on the small scales mesh
      dx = smallScales%mesh(2,1,1,1) - smallScales%mesh(1,1,1,1)
      dy = smallScales%mesh(1,2,1,2) - smallScales%mesh(1,1,1,2)
      dz = smallScales%mesh(1,1,2,3) - smallScales%mesh(1,1,1,3)
      
!TODO: This assumes cell data. Need to make this general for problems with edge
!data

      xst = smallScales%mesh(1,1,1,1) - dx/2
      yst = smallScales%mesh(1,1,1,2) - dy/2
      zst = smallScales%mesh(1,1,1,3) - dz/2

      xen = smallScales%mesh(smallScales%gpC%xsz(1),1,1,1) + dx/2
      yen = smallScales%mesh(1,smallScales%gpC%xsz(2),1,2) + dy/2
      zen = smallScales%mesh(1,1,smallScales%gpC%xsz(3),3) + dz/2

      ! Get the number of QH regions on the current MPI proc based on the local
      ! domain to domain size ratio
      this%isz = ceiling((xen - xst)/(xDom(2) - xDom(1))*this%nx)
      this%jsz = ceiling((yen - yst)/(yDom(2) - yDom(1))*this%ny)
      this%ksz = ceiling((zen - zst)/(zDom(2) - zDom(1))*this%nz)

      ! Now compute the global start and end indices
      call MPI_Allgather(this%isz,1,MPI_INTEGER,iszAll,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
      call MPI_Allgather(this%jsz,1,MPI_INTEGER,jszAll,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
      call MPI_Allgather(this%ksz,1,MPI_INTEGER,kszAll,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)

      this%ist = sum(iszAll(1:nrank)) + 1
      this%jst = sum(kszAll(1:nrank)) + 1
      this%kst = sum(jszAll(1:nrank)) + 1
    
      ! Check to see if nx, ny, or nz should be updated
      call assert(p_sum(this%isz*this%jsz*this%ksz) == ntotal,'Number of QH'//&
        ' regions has changed. This needs to be implemented')

      ! Allocate memory
      allocate(this%xC(this%isz), this%yC(this%jsz), this%zC(this%ksz))
      allocate(this%xE(this%isz+1), this%yE(this%jsz+1), this%zE(this%ksz+1))
      allocate(this%gID(this%isz,this%jsz,this%ksz))


      ! Define the mesh        
      this%xE = linspace(xst,xen,this%isz+1)
      this%yE = linspace(yst,yen,this%jsz+1) 
      this%zE = linspace(zst,zen,this%ksz+1)
!      this%xE = getMeshEdgeValues(this%ist,this%isz+1,this%dx,xDom(1))
!      this%yE = getMeshEdgeValues(this%jst,this%jsz+1,this%dy,yDom(1))
!      this%zE = getMeshEdgeValues(this%kst,this%ksz+1,this%dz,zDom(1))
      
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
    end subroutine
    
    subroutine destroy(this)
      class(QHmesh), intent(inout) :: this

      if (allocated(this%xE)) deallocate(this%xE)
      if (allocated(this%yE)) deallocate(this%yE)
      if (allocated(this%zE)) deallocate(this%zE)
      if (allocated(this%xC)) deallocate(this%xC)
      if (allocated(this%yC)) deallocate(this%yC)
      if (allocated(this%zC)) deallocate(this%zC)
    end subroutine
    
end module
