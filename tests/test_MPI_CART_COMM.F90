#include "../problems/gabor/hitEnrich_files/initialize.F90"

program test_MPI_CART_COMM
  use kind_parameters,    only: clen, rkind
  use mpi
  use IncompressibleGrid, only: igrid
  use enrichmentMod,      only: enrichmentOperator, testMPIsendRecv
  implicit none

  integer :: ierr
  character(len=clen) :: inputfileLS, inputfileSS, inputfileGM 
  type(igrid) :: largeScales, smallScales
  type(enrichmentOperator) :: enrich

  call MPI_Init(ierr)

  call GETARG(1,inputfileLS)
  call GETARG(2,inputfileSS)
  call GETARG(3,inputfileGM)

  call largeScales%init(inputfileLS, .true.) 
  call smallScales%init(inputfileSS, .false.)

  call enrich%init(smallScales,largeScales,inputfileGM,testOnly=.true.)

  call MPI_Finalize(ierr)
end program test_MPI_CART_COMM
