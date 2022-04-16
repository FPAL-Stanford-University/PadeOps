include "../problems/gabor/turbHalfChannel_files/initialize.F90"

program test_gaborMode_extrapolateToBoundary
  use kind_parameters, only: rkind, clen
  use largeScalesMod, only: getLargeScaleData, computeLargeScaleGradient
  use domainSetup, only: gpQHcent, getStartAndEndIndices
  use fortran_assert, only: assert
  use gaborIO_mod, only: readFields
  implicit none

  character(len=clen) :: inputfile, datadir, fname
  integer :: ierr, ioUnit
  real(rkind), dimension(:,:,:), allocatable :: L, KE
  integer :: ist, ien, jst, jen, kst, ken
  integer :: isz, jsz, ksz
  
  namelist /IO/ datadir

  ! Initialize MPI, HDF5, & generate numerical mesh
  call initializeProblem(inputfile)
 
  ! Get data directory path from inputfile
  open(unit=ioUnit, file=trim(inputfile), form='FORMATTED', iostat=ierr)
  read(unit=ioUnit, NML=IO)
  close(ioUnit)

  ! Read large scale data  
  write(fname,'(A)')trim(datadir)//'/'//'LargeScaleVelocityB.h5'
  call getLargeScaleData(trim(fname),.false.)

  ! Compute large scale gradient
  write(fname,'(A)')trim(datadir)//'/'//'LargeScaleGradientB.h5'
  call computeLargeScaleGradient(trim(fname))

  ! Read large scale parameters

    call getStartAndEndIndices(gpQHcent,ist,ien,jst,jen,kst,ken,isz,jsz,ksz)
    
    ! Allocate memory for integral scale and kinetic energy in each QH region
    allocate(L(isz,jsz,ksz), KE(isz,jsz,ksz))

    ! Read large scale parameters
    write(fname,'(A)')trim(datadir)//'/'//'LargeScaleParams.h5'
    call readFields(trim(fname),L,KE,'/L','/KE',gpQHcent)

    ! Free up memory
    deallocate(L,KE)
    call finalizeProblem()
end program
