program test_gaborMode_interpToMode
  use kind_parameters, only: rkind, clen
  use mpi
  use domainSetup, only: setupDomainXYperiodic, finalizeDomainSetup, &
    getStartAndEndIndices, gpLESe, xLESe, yLESe, zLESe, dxLES, dyLES, dzLES
  use GaborModeRoutines, only: gmxloc, gmyloc, gmzloc, interpToMode, &
    finalizeGaborModes
  use fortran_assert, only: assert
  use exits, only: message
  implicit none

  character(len=clen) :: inputfile
  integer :: ierr
  integer :: ist, ien, jst, jen, kst, ken, isz, jsz, ksz
  integer :: i, j, k
  real(rkind) :: tol = 1.d-13

  ! Data array defined on LES edge mesh
  real(rkind), dimension(:,:,:), allocatable :: datIn

  ! Interplated data
  real(rkind) :: datOut
  real(rkind) :: datCheck

  ! Interpolant coefficients
  real(rkind) :: a, b, c, d

  ! Initialize MPI
  call MPI_Init(ierr)

  ! Get inputfile path & name from command line  
  call GETARG(1,inputfile)

  ! Setup the domain
  call setupDomainXYperiodic(inputfile)

  ! Define data array
  a = 0.8147d2
  b = 0.9058d2
  c = 0.1270d2
  d = 0.9134d2

  call getStartAndEndIndices(gpLESe,ist,ien,jst,jen,kst,ken,isz,jsz,ksz)
  allocate(datIn(isz,jsz,ksz))
  do k = 1,ksz
    do j = 1,jsz
      do i = 1,isz
        datIn(i,j,k) = a*xLESe(i) + b*yLESe(j) + c*zLESe(k) + d
      end do
    end do
  end do

  ! Define mode locations
  allocate(gmxloc(1), gmyloc(1), gmzloc(1))
  gmxloc = 0.5469d0*(xLESe(isz)-xLESe(1)) + xLESe(1)
  gmyloc = 0.9575d0*(yLESe(jsz)-yLESe(1)) + yLESe(1)
  gmzloc = 0.9649d0*(zLESe(ksz)-zLESe(1)) + zLESe(1)

  call interpToMode(datIn,datOut,1,gpLESe,dxLES,dyLES,dzLES,xLESe,yLESe,zLESe)
  datCheck = a*gmxloc(1) + b*gmyloc(1) + c*gmzloc(1) + d

  call assert(abs(datOut - datCheck) < tol)
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  call message("Test PASSED!")

  ! Deallocate memory
  call finalizeGaborModes()
  call finalizeDomainSetup() 

  ! Finalize MPI
  call MPI_Finalize(ierr)
end program
