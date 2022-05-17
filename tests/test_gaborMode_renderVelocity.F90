#include "../problems/gabor/turbHalfChannel_files/turbHalfChanMod.F90"

program test_gaborMode_renderVelocity
  use kind_parameters, only: rkind, clen
  use largeScalesMod, only: readFields, writeFields
  use domainSetup, only: nxF, nyF, nzF, gpFE, getStartAndEndIndices, gpFC
  use GaborModeRoutines, only: renderVelocityXYperiodic, uhatR, uhatI, vhatR, vhatI, &
    whatR, whatI, kx, ky, kz, gmxloc, gmyloc, gmzloc, modeMemoryInitialized, &
    isotropicModesInitialized, uG, vG, wG, uGout, vGout, wGout
  use GaborIO_mod, only: readModes
  use fortran_assert, only: assert
  use mpi
  use decomp_2d, only: nrank
  use basic_io, only: read_1d_ascii
  use turbHalfChanMod, only: computeLargeScaleParams, initializeProblem, &
    finalizeProblem
  use exits, only: message
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
  !TODO: implement mode exchange using MPI
  !call modeExchange()
  isotropicModesInitialized = .true.

  ! Render the velocity field
  write(fname,'(A)')trim(datadir)//'/'//'IsotropicVelocity'
  if (newVelocity) then
    call renderVelocityXYperiodic()
    call writeFields(trim(fname)//'Fortran.h5',uG,vG,wG,'/u','/v','/w',gpFE)
    call writeFields(trim(fname)//'InterpFortran.h5',uGout,vGout,wGout,&
      '/u','/v','/w',gpFC)
  else
    call readFields(trim(fname)//'Fortran.h5',uG,vG,wG,'/u','/v','/w',gpFE)
    call readFields(trim(fname)//'FortranInterp.h5',uGout,vGout,wGout,&
      '/u','/v','/w',gpFE)
  end if

  ! Read comparison velocity data
  call getStartAndEndIndices(gpFE,ist,ien,jst,jen,kst,ken,isz,jsz,ksz)
  allocate(uTrue(ist:ien,jst:jen,kst:ken))
  allocate(vTrue(ist:ien,jst:jen,kst:ken))
  allocate(wTrue(ist:ien,jst:jen,kst:ken))
  write(fname,'(A)')trim(datadir)//'/'//'IsotropicVelocityMATLAB.h5'
  call readFields(trim(fname),uTrue,vTrue,wTrue,'/u','/v','/w',gpFE)
  
  call assert(maxval(abs(uTrue - uG)) < small,'u discrepancy')
  call assert(maxval(abs(vTrue - vG)) < small,'v discrepancy')
  call assert(maxval(abs(wTrue - wG)) < small,'w discrepancy')

  deallocate(uTrue,vTrue,wTrue)
  call getStartAndEndIndices(gpFC,ist,ien,jst,jen,kst,ken,isz,jsz,ksz)
  allocate(uTrue(ist:ien,jst:jen,kst:ken))
  allocate(vTrue(ist:ien,jst:jen,kst:ken))
  allocate(wTrue(ist:ien,jst:jen,kst:ken))
  write(fname,'(A)')trim(datadir)//'/'//'IsotropicVelocityInterpMATLAB.h5'
  call readFields(trim(fname),uTrue,vTrue,wTrue,'/u','/v','/w',gpFC)

  call message('Comparing interpolated fields:') 
  print*, "maxval(abs(uTrue - uGout)):", maxval(abs(uTrue - uGout)) 
  print*, "maxval(abs(vTrue - vGout)):", maxval(abs(vTrue - vGout)) 
  print*, "maxval(abs(wTrue - wGout)):", maxval(abs(wTrue - wGout)) 
  
  !call assert(maxval(abs(uTrue - uGout)) < small,'u discrepancy')
  !call assert(maxval(abs(vTrue - vGout)) < small,'v discrepancy')
  !call assert(maxval(abs(wTrue - wGout)) < small,'w discrepancy')

  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  if (nrank == 0) print*, "Test PASSED!"
  call MPI_Barrier(MPI_COMM_WORLD,ierr)

  ! Finalize problem
  deallocate(uTrue,vTrue,wTrue)
  call finalizeProblem()
end program
