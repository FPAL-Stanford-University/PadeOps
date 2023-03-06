#include "../problems/incompressible/randomJetHIT_files/jetMod.F90"
program test_randomJetHITrandomSampling
  use kind_parameters, only: rkind, clen
  use jetMod,          only: jet, outputdir
  use mpi
  use constants,       only: pi
  use basic_io,        only: write_1d_ascii, write_2d_ascii
  use timer,           only: tic, toc
  implicit none

  type(jet), dimension(:,:), allocatable :: jetArray
  character(len=clen) :: inputfile, fname
  integer :: Nsteps = 1000000
  integer :: ierr, i, j, gID, NjetsTotal, n, realization, iter
  real(rkind) :: Lx = 4.d0*pi, Ly = 4.d0*pi, Lz
  integer :: NjetsX = 8, NjetsY = 8
  real(rkind) :: jetXsz = 0.5d0, jetYsz = 0.5d0
  real(rkind) :: dxJet, dyJet, xloc, yloc, time
  real(rkind), parameter :: dt = 1.d-3
  logical :: myStatus
  real(rkind) :: jetZsz, jetZst
  integer :: ntimesFilter, dumpMaskFreq, ioUnit
  real(rkind), dimension(:), allocatable :: jetID, duration
  real(rkind), dimension(:,:), allocatable :: location
  integer :: dumpFreq = 100

  call MPI_Init(ierr)
  call getarg(1,inputfile)
  
  namelist /PROBLEM_INPUT/ Lx, Ly, Lz, ntimesFilter, dumpMaskFreq, &
    NjetsX, NjetsY, jetZst, jetXsz, jetYsz, jetZsz, dumpFreq, Nsteps

  ioUnit = 11
  open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
  read(unit=ioUnit, NML=PROBLEM_INPUT)
  close(ioUnit)    
  
  allocate(jetArray(NjetsX,NjetsY))
  
  dxJet = Lx/real(NjetsX,rkind)
  dyJet = Ly/real(NjetsY,rkind)
  
  NjetsTotal = NjetsY*NjetsX

  do j = 1,NjetsY
    yloc = 0.5d0*dyJet + real(j - 1,rkind)*dyJet
    do i = 1,NjetsX
      xloc = 0.5d0*dxJet + real(i - 1,rkind)*dxJet
      gID = (j-1)*NjetsX + i
      call jetArray(i,j)%init(xloc,yloc,0.d0,jetXsz,jetYsz,2.d0,gID,&
        inputfile,NjetsTotal,0.d0)
    end do
  end do

  allocate(jetID(NjetsX*NjetsY))
  allocate(duration(NjetsX*NjetsY))
  allocate(location(NjetsX*NjetsY,2))
  realization = 0
  call tic()
  do n = 1,Nsteps
    time = real(n-1)*dt
    iter = 0
    do j = 1, NjetsY
      do i = 1,NjetsX
        myStatus = jetArray(i,j)%isOn(time)
        if (mod(n,dumpFreq) == 0 .and. myStatus) then
          iter = iter + 1
          jetID(iter)      = real(jetArray(i,j)%getID())
          duration(iter)   = jetArray(i,j)%getDuration()
          location(iter,:) = jetArray(i,j)%getLocation()
        end if
      end do
    end do
    if (mod(n,dumpFreq) == 0) then
      call toc(n,Nsteps,whichMssg = 2)
      realization = realization + 1
      write(fname,'(A,I6.6,A)')trim(outputdir)//'jetID_',realization,'.txt'
      call write_1d_ascii(jetID(1:iter),trim(fname))
      write(fname,'(A,I6.6,A)')trim(outputdir)//'location_',realization,'.txt'
      call write_2d_ascii(location(1:iter,:),trim(fname))
      write(fname,'(A,I6.6,A)')trim(outputdir)//'duration_',realization,'.txt'
      call write_1d_ascii(duration(1:iter),trim(fname))
    end if
  end do

  do j = 1, NjetsY
    do i = 1,NjetsX
      call jetArray(i,j)%destroy()
    end do
  end do
  
  deallocate(jetID, location, duration)
  call MPI_finalize(ierr)
end program
