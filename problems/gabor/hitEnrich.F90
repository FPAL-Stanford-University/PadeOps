#include "hitEnrich_files/initialize.F90"

program hitEnrich
  use kind_parameters,         only: clen
  use enrichmentMod,           only: enrichmentOperator
  use IncompressibleGrid,      only: igrid
  use auxiliary_openmp_subs,   only: GetArguments
  use mpi
  use reductions
  use decomp_2d,               only: nrank, nproc
  implicit none

  character(len=clen) :: inputfileLS, inputfileSS, inputfileGM 
  integer :: ierr, provided
  type(igrid) :: largeScales, smallScales
  type(enrichmentOperator) :: enrich
  integer :: n
  integer :: tid = 0
      
  call MPI_Init_thread(MPI_THREAD_FUNNELED,provided,ierr)
  
  call GETARG(1,inputfileLS)
  call GETARG(2,inputfileSS)
  call GETARG(3,inputfileGM)
  
  call largeScales%init(inputfileLS, .true.) 
  call smallScales%init(inputfileSS, .false.)

  call enrich%init(smallScales,largeScales,inputfileGM)

  call largeScales%initLargeScales(tid,readGradients=.false.)

  call enrich%generateModes()
  call enrich%renderVelocity()
  call enrich%dumpSmallScales()
  print*, p_sum(sum(abs(smallScales%u))) 
  
  !do while (enrich%continueSimulation())
  !  call enrich%updateLargeScales(timeAdvance=.true.)
  !  call enrich%advanceTime()
  !  call enrich%wrapupTimeStep()
  !end do

  call enrich%destroy()
  call largeScales%destroy()
  call smallScales%destroy()

  call MPI_Finalize(ierr)
end program hitEnrich
