#include "turbHalfChannel_files/turbHalfChanMod.F90"

program turbHalfChannel
  use turbHalfChanMod, only: initializeProblem, finalizeProblem, &
    computeLargeScaleParams
  use largeScalesMod, only: getLargeScaleData, computeLargeScaleGradient, &
    writeFields
  use gaborModeRoutines, only: generateIsotropicModes, renderVelocityXYperiodic,&
    uGout, vGout, wGout, strainIsotropicModes
  use domainSetup, only: gpFC
  use kind_parameters, only: clen
  
  implicit none
  character(len=clen) :: inputfile, datadir, fname, outputdir
  integer :: ioUnit, ierr
        
  namelist /IO/ datadir, outputdir
  
  call initializeProblem(inputfile)

  ! Read inputfile
  ioUnit = 1
  open(unit=ioUnit, file=trim(inputfile), form='FORMATTED', iostat=ierr)
  read(unit=ioUnit, NML=IO)
  close(ioUnit)
  
  ! Read large scale initial condition data
  write(fname,'(A)')trim(datadir)//'/'//'LargeScaleField.h5'
  call getLargeScaleData(fname)

  ! TODO: Extend large scale data to domain boundaries

  ! Compute or read large scale parameters
  write(fname,'(A)')trim(datadir)//'/'//'LargeScaleParams.h5'
  call computeLargeScaleParams(fname,'/KE','/L') 

  ! TODO: extend large scale fields to domain boundaries

  ! Initialize the Gabor Modes
  call generateIsotropicModes()
  call strainIsotropicModes()

  ! Render initialized velocity field (optional)
  call renderVelocityXYperiodic()
  write(fname,'(A)')trim(outputdir)//'/'//'GaborInducedVelocity.h5'
  call writeFields(fname,uGout,vGout,wGout,'/u','/v','/w',gpFC)

  call finalizeProblem()
end program turbHalfChannel