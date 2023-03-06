#include "../problems/incompressible/randomJetHIT_files/jetMod.F90"
program test_randomJetHITrandomSampling
  use kind_parameters, only: rkind, clen
  use jetMod,          only: jet
  use mpi
  use constants,       only: pi
  implicit none

  type(jet), dimension(:,:), allocatable :: jetArray
  character(len=clen) :: inputfile
  integer, parameter :: Nsamples = 10000
  integer :: ierr, i, j, gID, NjetsTotal
  real(rkind), parameter :: Lx = 4.d0*pi, Ly = 4.d0*pi
  integer, parameter :: NjetsX = 8, NjetsY = 8
  real(rkind) :: dxJet, dyJet, xloc, yloc

  call MPI_Init(ierr)
  call getarg(1,inputfile)
  
  allocate(jetArray(NjetsX,NjetsY))
  
  dxJet = Lx/real(NjetsX,rkind)
  dyJet = Ly/real(NjetsY,rkind)
  
  NjetsTotal = NjetsY*NjetsX

  do j = 1,NjetsY
    yloc = 0.5d0*dyJet + real(j - 1,rkind)*dyJet
    do i = 1,NjetsX
      xloc = 0.5d0*dxJet + real(i - 1,rkind)*dxJet
      gID = (j-1)*NjetsX + i
      call jetArray(i,j)%init(xloc,yloc,0.d0,0.5d0,0.5d0,2.d0,gID,&
        inputfile,NjetsTotal,0.d0)
    end do
  end do

  do j = 1, NjetsY
    do i = 1,NjetsX
      call jetArray(i,j)%testRandomSampling(Nsamples)
    end do
  end do
  do j = 1, NjetsY
    do i = 1,NjetsX
      call jetArray(i,j)%destroy()
    end do
  end do

  call MPI_finalize(ierr)
end program
