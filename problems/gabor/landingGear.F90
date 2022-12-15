#include "landingGear_files/initialize.F90"

program landingGear
  use kind_parameters,         only: clen
  use enrichmentMod,           only: enrichmentOperator, nthreads
  use IncompressibleGrid,      only: igrid
  use auxiliary_openmp_subs,   only: GetArguments
  use mpi
  implicit none

  character(len=clen) :: inputS1, inputS2, inputS3, inputG12, inputG23
  integer :: ierr, provided
  type(igrid) :: S1, S2
  type(enrichmentOperator) :: enrich12 
      
  call MPI_Init_thread(MPI_THREAD_FUNNELED,provided,ierr)
  
  call GETARG(1,inputS1)
  call GETARG(2,inputS2)
  call GETARG(3,inputG12)
  call GetArguments(nthreads,4)

  call S1%init(inputS1, .true.) 
  call S2%init(inputS2, .false.)
  
  call enrich12%init(S2,S1,inputG12)
  call enrich12%renderVelocity()
  call enrich12%dumpSmallScales()

  call enrich12%destroy()
  call S1%destroy()
  call S2%destroy()

  call MPI_Finalize(ierr)
end program landingGear
