program test_gaborMode_ffts
  use mpi
  use decomp_2d
  use fftstuff,        only: ffts
  use fortran_assert,  only: assert
  use kind_parameters, only: rkind
  use constants,       only: pi
  implicit none

  type(ffts) :: fftx, ffty, fftz
  integer :: ierr
  integer, parameter :: nx = 64, ny = 64, nz = 64
  real(rkind) :: Lx = 2.d0*pi, Ly = 2.d0*pi, Lz = 2.d0*pi
  real(rkind), dimension(:,:,:), allocatable :: x, y, z, f, df
  complex(rkind), dimension(:,:,:), allocatable :: fhat

  ! Initialize things
  call MPI_Init(ierr)
  call decomp_2d_init(nx,ny,nz,0,0)
  
  ! Allocate arrays
  allocate(f(xsize(1),xsize(2),xsize(3)))
  allocate(df(xsize(1),xsize(2),xsize(3)))
  allocate(x(xsize(1),xsize(2),xsize(3)))
  allocate(y(xsize(1),xsize(2),xsize(3)))
  allocate(z(xsize(1),xsize(2),xsize(3)))

  ! Grid spacing
  dx = Lx/real(nx,rkind)
  dy = Ly/real(ny,rkind)
  dz = Lz/real(nz,rkind)

  ! Mesh gen
  do k = 1,nz
    do j = 1,ny
      do i = 1,nx
        x(i,j,k) = (i-1)*dx
        y(i,j,k) = (j-1)*dy
        z(i,j,k) = (k-1)*dz
      end do
    end do
  end do

  ! Initialize ffts
  ierr =  fftx%init(xsize(1),'x',xsize(2),xsize(3),dx)
  call assert(ierr == 0,'fftx init')
  ierr =  ffty%init(ysize(2),'y',ysize(1),ysize(3),dy)
  call assert(ierr == 0,'ffty init')
  ierr =  fftz%init(zsize(3),'z',zsize(1),zsize(2),dz)
  call assert(ierr == 0,'fftz init')

  ! Allocate complex array
  allocate(fhat(fftx%split,fftx%n2,fftx%n3))

  ! Define input function
  f = sin(2.d0*x)*sin(2.d0*y)*sin(2.d0*z)

  ! Compute derivative in spectral space

  ! Invert the transform

  ! Compare to analytical derivative

  ! Finalize things
  call decomp_2d_finalize()
  call MPI_Finalize(ierr)
end program
