#include "../problems/gabor/turbHalfChannel_files/turbHalfChanMod.F90"

program test_gaborMode_IsotropicModes
  use kind_parameters, only: rkind, clen
  use largeScalesMod, only: getLargeScaleData, computeLargeScaleGradient, U
  use domainSetup, only: yLES
  use GaborModeRoutines, only: generateIsotropicModes, getModelSpectrum, nk, &
    renderVelocity
  use GaborIO_mod, only: writeModes, readModes
  use fortran_assert, only: assert
  use mpi
  use decomp_2d, only: nrank, nproc
  use gridtools, only: logspace
  use basic_io, only: read_1d_ascii
  use turbHalfChanMod, only: computeLargeScaleParams, initializeProblem, &
    finalizeProblem
  implicit none
  character(len=clen) :: inputfile, fname, datadir, outputdir
  integer :: ioUnit, ierr
  real(rkind), dimension(:), allocatable :: Etrue, E, k
  real(rkind) :: small = 1.d-12

  namelist /IO/ datadir, outputdir

  ! Initiaize MPI, HDF5, & generate numerical mesh
  call initializeProblem(inputfile)

  ! Read inputfile to get datadir
  ioUnit = 1
  open(unit=ioUnit, file=trim(inputfile), form='FORMATTED', iostat=ierr)
  read(unit=ioUnit, NML=IO)
  close(ioUnit)
 
  ! Get large scale data 
  write(fname,'(A)')trim(datadir)//'/'//'LargeScaleVelocity.h5'
  call getLargeScaleData(trim(fname))

  write(fname,'(A)')trim(datadir)//'/'//'LargeScaleGradient.h5'
  call computeLargeScaleGradient(trim(fname))

  write(fname,'(A)')trim(datadir)//'/'//'LargeScaleParams.h5'
  call computeLargeScaleParams(trim(fname),'/KE','/L')

  ! Test model spectrum function
  allocate(E(nk), k(nk))
  k = logspace(log10(8.d0),log10(32.d0),nk)
  E = getModelSpectrum(k,4.2d0,3.221d0,nk)
  write(fname,'(A)')trim(datadir)//'/'//'Emodel.dat'
  call read_1d_ascii(Etrue,trim(fname))
  call assert(maxval(abs(E-Etrue)) < small,'E discrepency')
  deallocate(E,k,Etrue)
  
  ! Initialize isotropic Gabor modes
  call generateIsotropicModes()
  write(fname,'(A,I1,A)')trim(outputdir)//'/'//'GaborModes',nproc,'PE.h5'
  call writeModes(trim(fname))

  ! Render the velocity field
  call renderVelocity()

  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  if (nrank == 0) print*, "Test PASSED!"
  call MPI_Barrier(MPI_COMM_WORLD,ierr)

  ! Finalize problem
  call finalizeProblem()
end program
