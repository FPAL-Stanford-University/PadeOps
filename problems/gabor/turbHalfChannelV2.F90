#include "turbHalfChannel_files/turbHalfChanMod.F90"

program turbHalfChannel
  !use turbHalfChanMod, only: initializeProblem, finalizeProblem, &
  !  computeLargeScaleParams
  !use largeScalesMod, only: getLargeScaleData, computeLargeScaleGradient, &
  !  writeFields
  !use gaborModeRoutines, only: generateIsotropicModes, renderVelocityXYperiodic,&
  !  uGout, vGout, wGout, strainIsotropicModes
  !use domainSetup, only: gpFC
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
  
  call largeScales%init(inputfileLS, .true.) 
  call smallScales%init(inputfileSS, .false.)

  call enrich%init(smallScales,largeScales,inputfileGM)

  !do tid = 1,nsteps
  !  call enrichmentOperator%updateLargeScales()
  !  call enrichmentOperator%advanceTime(dt,...)
  !  if (mod(tid,trender) == 0) then
  !    call enrichmentOperator%render()
  !  end if
  !  if (mod(tid,twrite) == 0) then
  !    call smallScales%dump_data()
  !  end if
  !end do

  call enrich%destroy()
  call largeScales%destroy()
  call smallScales%destroy()

  
  call MPI_Finalize(ierr)
end program turbHalfChannel
