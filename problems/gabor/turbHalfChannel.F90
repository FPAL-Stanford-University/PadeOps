include 'turbulenHalfChannel_files/intialize.F90'

program turbHalfChannel
  use largeScalesMod, only: getLargeScaleData, computeLargeScaleGradient
  
  implicit none
  character(len=clen) :: inputfile
  integer :: ierr
  
  call initializeProblem(inputfile)

  ! Read large scale initial condition data
  write(fname,'(A)')trim(datadir)//'/'//'largeScaleFile.h5'
  call getLargeScaleData(fname)

  ! Compute large scale gradient from velocity data
  call computeLargeScaleGradient()

  ! TODO: Extend large scale data to domain boundaries

  ! TODO: Compute or read large scale parameters

  ! TODO: extend large scale fields to domain boundaries

  call finalizeProblem()
end program turbHalfChannel
