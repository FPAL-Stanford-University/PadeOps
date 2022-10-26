#include "hitEnrich_files/initialize.F90"

program hitEnrich
  use kind_parameters,         only: clen
  use enrichmentMod,           only: enrichmentOperator, nthreads
  use IncompressibleGrid,      only: igrid
  use HIT_Periodic_parameters, only: Lx, Ly, Lz
  use auxiliary_openmp_subs,   only: GetArguments
  use mpi
  implicit none

  character(len=clen) :: inputfileLS, inputfileSS, inputfileGM 
  character(len=clen) :: datadir, fname, outputdir
  integer :: ioUnit, ierr, provided
  type(igrid) :: largeScales, smallScales
  type(enrichmentOperator) :: enrich
      
  call MPI_Init_thread(MPI_THREAD_FUNNELED,provided,ierr)
  
  call GETARG(1,inputfileLS)
  call GETARG(2,inputfileSS)
  call GETARG(3,inputfileGM)
  call GetArguments(nthreads)
  
  call largeScales%init(inputfileLS, .true.) 
  call smallScales%init(inputfileSS, .false.)

  call enrich%init(smallScales,largeScales,inputfileGM,Lx,Ly,Lz)

  do while (enrich%continueSimulation())
    call enrich%updateLargeScales(timeAdvance=.true.)
    call enrich%advanceTime()
    call enrich%wrapupTimeStep()
  end do

  call enrich%destroy()
  call largeScales%destroy()
  call smallScales%destroy()

  call MPI_Finalize(ierr)
end program hitEnrich
