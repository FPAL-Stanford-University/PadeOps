subroutine initializeProblem(inputfile)
  use domainSetup, only: setupDomainXYperiodic
  use largeScalesMod, only: initLargeScales
  use kind_parameters, only: clen
  implicit none
  character(len=clen), intent(inout) :: inputfile
  call initializeExternalLibraries()
  call GETARG(1,inputfile)
  call setupDomainXYperiodic(inputfile)
  call initLargeScales()
end subroutine

subroutine finalizeProblem()
  use domainSetup, only: finalizeDomainSetup
  use largeScalesMod, only: finalizeLargeScales
  implicit none
  call finalizeDomainSetup()
  call finalizeLargeScales()
  call finalizeExternalLibraries()
end subroutine

subroutine initializeExternalLibraries()
  use fortran_assert, only: assert
  use mpi
  use hdf5
  implicit none
  integer :: ierr

  call MPI_Init(ierr)
  call assert(ierr == MPI_SUCCESS,'MPI initialization failed.')
  call H5open_f(ierr)
  call assert(ierr == 0,'HDF5 initialization failed.')
end subroutine

subroutine finalizeExternalLibraries()
  use fortran_assert, only: assert
  use mpi
  use hdf5
  implicit none
  integer :: ierr

  call H5close_f(ierr)
  call assert(ierr == 0,'HDF5 finalization failed.')
  call MPI_Finalize(ierr)
  call assert(ierr == MPI_SUCCESS,'MPI finalization failed.')
end subroutine
