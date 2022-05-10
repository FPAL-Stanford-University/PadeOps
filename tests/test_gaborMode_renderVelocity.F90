#include "../problems/gabor/turbHalfChannel_files/turbHalfChanMod.F90"

program test_gaborMode_renderVelocity
  use kind_parameters, only: rkind, clen
  use largeScalesMod, only: readFields, writeFields
  use domainSetup, only: nxF, nyF, nzF, gpFE, getStartAndEndIndices
  use GaborModeRoutines, only: renderVelocityXYperiodic, uhatR, uhatI, vhatR, vhatI, &
    whatR, whatI, kx, ky, kz, gmxloc, gmyloc, gmzloc, modeMemoryInitialized, &
    isotropicModesInitialized, uG, vG, wG
  use GaborIO_mod, only: readModes
  use fortran_assert, only: assert
  use mpi
  use decomp_2d, only: nrank
  use basic_io, only: read_1d_ascii
  use turbHalfChanMod, only: computeLargeScaleParams, initializeProblem, &
    finalizeProblem
  implicit none
  character(len=clen) :: inputfile, fname, datadir, outputdir
  integer :: ioUnit, ierr
  real(rkind), dimension(:,:,:), allocatable :: uTrue, vTrue, wTrue
  integer :: ist, ien, jst, jen, kst, ken
  integer :: isz, jsz, ksz
  real(rkind) :: small = 1.d-13
  logical :: newVelocity = .true.

  namelist /IO/ datadir, outputdir

  ! Initiaize MPI, HDF5, & generate numerical mesh
  call initializeProblem(inputfile)

  ! Read inputfile to get datadir
  ioUnit = 1
  open(unit=ioUnit, file=trim(inputfile), form='FORMATTED', iostat=ierr)
  read(unit=ioUnit, NML=IO)
  close(ioUnit)
 
  ! Read Gabor mode info
  call assert(modeMemoryInitialized,"Gabor mode memory not allocated - TEST")
  write(fname,'(A)')trim(datadir)//'/'//'IsotropicModesMATLAB.h5'
  call readModes(fname,uhatR,uhatI,vhatR,vhatI,whatR,whatI,kx,ky,kz,gmxloc,&
    gmyloc,gmzloc)
  isotropicModesInitialized = .true.

  ! Render the velocity field
  write(fname,'(A)')trim(datadir)//'/'//'IsotropicVelocityFortran.h5'
  if (newVelocity) then
    call renderVelocityXYperiodic()
    call writeFields(trim(fname),uG,vG,wG,'/u','/v','/w',gpFE)
  else
    call readFields(trim(fname),uG,vG,wG,'/u','/v','/w',gpFE)
  end if

  ! Read comparison velocity data
  call getStartAndEndIndices(gpFE,ist,ien,jst,jen,kst,ken,isz,jsz,ksz)
  !call getStartAndEndIndices(gpFh,ist,ien,jst,jen,kst,ken,isz,jsz,ksz)
  allocate(uTrue(ist:ien,jst:jen,kst:ken))
  allocate(vTrue(ist:ien,jst:jen,kst:ken))
  allocate(wTrue(ist:ien,jst:jen,kst:ken))
  write(fname,'(A)')trim(datadir)//'/'//'IsotropicVelocityMATLAB.h5'
  call readFields(trim(fname),uTrue,vTrue,wTrue,'/u','/v','/w',gpFE)
  !write(fname,'(A)')trim(datadir)//'/'//'IsotropicVelocityHaloMATLAB.h5'
  !call readFields(trim(fname),uTrue,vTrue,wTrue,'/u','/v','/w',gpFh)
  
  print*, "maxval(abs(uTrue)):", maxval(abs(uTrue))
  print*, "maxval(abs(uG)):", maxval(abs(uG))
  call assert(maxval(abs(uTrue - uG)) < small,'u discrepancy')
  call assert(maxval(abs(vTrue - vG)) < small,'v discrepancy')
  call assert(maxval(abs(wTrue - wG)) < small,'w discrepancy')
  !print*, "maxval(abs(uFh)):", maxval(abs(uFh))
  !call assert(maxval(abs(uTrue - uFh)) < small,'u discrepancy')
  !call assert(maxval(abs(vTrue - vFh)) < small,'v discrepancy')
  !call assert(maxval(abs(wTrue - wFh)) < small,'w discrepancy')

  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  if (nrank == 0) print*, "Test PASSED!"
  call MPI_Barrier(MPI_COMM_WORLD,ierr)

  ! Finalize problem
  deallocate(uTrue,vTrue,wTrue)
  call finalizeProblem()
end program
