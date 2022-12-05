#include "hitEnrich_files/initialize.F90"

program hitEnrich
  use kind_parameters,         only: clen
  use enrichmentMod,           only: enrichmentOperator, nthreads
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
      
  call MPI_Init_thread(MPI_THREAD_FUNNELED,provided,ierr)
  
  call GETARG(1,inputfileLS)
  call GETARG(2,inputfileSS)
  call GETARG(3,inputfileGM)
  call GetArguments(nthreads,4)
  
  call largeScales%init(inputfileLS, .true.) 
  call smallScales%init(inputfileSS, .false.)

  call enrich%init(smallScales,largeScales,inputfileGM)

  print*, p_maxval(maxval(smallScales%u)) - p_minval(minval(smallScales%u))
! DEBUG
!do n = 1,nproc
!  if (nrank == n-1) then
!    print*, "PE", nrank
!    print*, enrich%QHgrid%gpC%xsz(1), enrich%QHgrid%gpC%xsz(2), enrich%QHgrid%gpC%xsz(3)
!    print*, enrich%QHgrid%gpC%xst(1), enrich%QHgrid%gpC%xst(2), enrich%QHgrid%gpC%xst(3)
!  end if
!  call MPI_Barrier(MPI_COMM_WORLD,ierr)
!end do
! END DEBUG
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
