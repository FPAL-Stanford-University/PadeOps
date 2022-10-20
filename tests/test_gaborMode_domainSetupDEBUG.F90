program domainSetupDEBUG
  use kind_parameters, only: rkind, clen
  use mpi
  use decomp_2d, only: nrank, nproc
  use domainSetup
  use fortran_assert, only: assert
  use exits, only: message
  implicit none

  integer :: ierr
  character(len=clen) :: inputfile, mssg
  integer :: ist1, ien1, jst1, jen1, kst1, ken1, isz1, jsz1, ksz1
  integer :: ist2, ien2, jst2, jen2, kst2, ken2, isz2, jsz2, ksz2
  real(rkind) :: tol = 1.d-13

  call MPI_Init(ierr)
  
  ! Get the input file path and file name
  call GETARG(1,inputfile)

  ! Setup the domain
  call setupDomainXYperiodic(inputfile)

  call getStartAndEndIndices(gpQHcent,ist1,ien1,jst1,jen1,kst1,ken1,isz1,jsz1,ksz1)
  call getStartAndEndIndices(gpLESe,ist2,ien2,jst2,jen2,kst2,ken2,isz2,jsz2,ksz2)
  write(mssg,'(A,F15.12,A,F15.12)')"yLESe(jsz2) = ", yLESe(jsz2), &
    " | yQHedge(jsz1+1) = ", yQHedge(jsz1+1)
  call assert(abs(yLESe(jsz2) - yQHedge(jsz1+1)) < tol,trim(mssg),nrank)

  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  call message('Test PASSED!')
  call MPI_Finalize(ierr)
end program
