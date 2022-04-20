#include 'turbulenHalfChannel_files/intialize.F90'
#include 'turbulenHalfChannel_files/computLargeScaleParams.F90'

program turbHalfChannel
  use largeScalesMod, only: getLargeScaleData, computeLargeScaleGradient, KE, L
  
  implicit none
  character(len=clen) :: inputfile, datadir
  integer :: ierr
  
  call initializeProblem(inputfile,datadir)

  ! Read large scale initial condition data
  write(fname,'(A)')trim(datadir)//'/'//'largeScaleField.h5'
  call getLargeScaleData(fname)

  ! Compute or read large scale gradient from velocity data
  call computeLargeScaleGradient()

  ! TODO: Extend large scale data to domain boundaries

  ! TODO: Compute or read large scale parameters
  write(fname,'()')trim(datadir)//'/'//'LargeScaleParams.h5'
  call computeLargeScaleParams(fname,KE,L,'/KE','/L') 

  ! TODO: extend large scale fields to domain boundaries

  call finalizeProblem()
end program turbHalfChannel
