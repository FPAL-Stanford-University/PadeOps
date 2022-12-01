#include "hitMultiScale_files/initialize.F90"

program hitMultiScale
  use kind_parameters,         only: clen
  use enrichmentMod,           only: enrichmentOperator, nthreads
  use IncompressibleGrid,      only: igrid
  use auxiliary_openmp_subs,   only: GetArguments
  use mpi
  implicit none

  character(len=clen) :: inputS1, inputS2, inputS3, inputG12, inputG23
  integer :: ierr, provided
  type(igrid) :: S1, S2, S3
  type(enrichmentOperator) :: enrich12, enrich23 
      
  call MPI_Init_thread(MPI_THREAD_FUNNELED,provided,ierr)
  
  call GETARG(1,inputS1)
  call GETARG(2,inputS2)
  call GETARG(3,inputS3)
  call GETARG(3,inputG12)
  call GETARG(3,inputG23)
  call GetArguments(nthreads)
  
  call S1%init(inputS1, .true.) 
  call S2%init(inputS2, .false.)
  call S3%init(inputS3, .false.)

  call enrich12%init(S1,S2,inputG12)
  call enrich23%init(S2,S3,inputG23)


  call enrich12%destroy()
  call enrich23%destroy()
  call S1%destroy()
  call S2%destroy()
  call S3%destroy()

  call MPI_Finalize(ierr)
end program hitMultiScale
