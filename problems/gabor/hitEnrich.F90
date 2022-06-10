#include "hitEnrich_files/initialize.F90"

program hitEnrich
  use kind_parameters, only: clen
  use enrichmentMod, only: enrichmentOperator  
  use IncompressibleGrid, only: igrid
  implicit none
  character(len=clen) :: inputfileLS, inputfileSS, inputfileGM 
  character(len=clen) :: datadir, fname, outputdir
  integer :: ioUnit, ierr
  type(igrid) :: largeScales, smallScales
  type(enrichmentOperator) :: enrich
      
  call MPI_Init(ierr)                                                
  
  call GETARG(1,inputfileLS)
  call GETARG(2,inputfileSS)
  call GETARG(3,inputfileGM)
  
  call largeScales%init(inputfileLS, .true.) 
  call smallScales%init(inputfileSS, .false.)

  call enrich%init(smallScales,largeScales,inputfileGM)

  do while (enrich%continueSimulation())
    call enrich%updateLargeScales()
    call enrich%advanceTime()
    call enrich%wrapupTimeStep()
  end do

  call enrich%destroy()
  call largeScales%destroy()
  call smallScales%destroy()

  
  call MPI_Finalize(ierr)
end program hitEnrich
