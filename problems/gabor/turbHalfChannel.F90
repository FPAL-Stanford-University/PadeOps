! Template for PadeOps

program turbHalfChannel
  use mpi
  
  implicit none

  integer :: ierr

  call MPI_Init(ierr)       ! Initialize MPI
  call GETARG(1,inputfile)  ! Get the location of the input file
  call setupDomainXYperiodic(inputfile) ! Define all numerical meshes (LES, QH, and high resolution)


  call MPI_Finalize(ierr) 
end program turbHalfChannel
