#include "../problems/incompressible/randomJetHIT_files/jetMod.F90"
program test_randomJetHITupdateWeight
  use kind_parameters, only: rkind, clen
  use jetMod,          only: jet
  use mpi
  use constants,       only: pi
  implicit none

  type(jet) :: J
  character(len=clen) :: inputfile
  integer :: ierr

  call MPI_Init(ierr)
  call getarg(1,inputfile)
  
  call J%init(1.d0,1.d0,0.d0,0.5d0,0.5d0,2.d0,1,inputfile,1,0.d0)
  call J%testUpdateWeight()
  call J%destroy()


  call MPI_finalize(ierr)
end program
