#include "../problems/gabor/turbHalfChannel_files/turbHalfChanMod.F90"

program test_gaborMode_strainModes
  use kind_parameters, only: rkind, clen
  use mpi
  use gaborModeRoutines, only: uhatR, uhatI, vhatR, vhatI, whatR, whatI, &
    kx, ky, kz, gmxloc, gmyloc, gmzloc, strainIsotropicModes 
  use turbHalfChanMod, only: initializeProblem, finalizeProblem, computeLargeScaleParams
  use largeScalesMod, only: getLargeScaleData 
  use fortran_assert, only: assert
  use basic_io, only: read_1d_ascii
  use exits, only: message
  implicit none

  integer :: ierr, ioUnit
  character(len=clen) :: inputfile, datadir, fname, outputdir
  real(rkind), dimension(:), allocatable :: uR, uI, vR, vI, wR, wI, k1, k2, k3, &
    xloc, yloc, zloc
  real(rkind) :: tol = 1.d-13
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

  ! Compute or read large scale parameters
  write(fname,'(A)')trim(datadir)//'/'//'LargeScaleParams.h5'
  call computeLargeScaleParams(fname,'/KE','/L') 

  call assert(size(gmxloc) == 25600)

  ! Read isotropic mode data
  write(fname,'(A)')trim(datadir)//'/'//'isotropicModes/uhatR.dat'
  call read_1d_ascii(uR,trim(fname))
  write(fname,'(A)')trim(datadir)//'/'//'isotropicModes/uhatI.dat'
  call read_1d_ascii(uI,trim(fname))
  write(fname,'(A)')trim(datadir)//'/'//'isotropicModes/vhatR.dat'
  call read_1d_ascii(vR,trim(fname))
  write(fname,'(A)')trim(datadir)//'/'//'isotropicModes/vhatI.dat'
  call read_1d_ascii(vI,trim(fname))
  write(fname,'(A)')trim(datadir)//'/'//'isotropicModes/whatR.dat'
  call read_1d_ascii(wR,trim(fname))
  write(fname,'(A)')trim(datadir)//'/'//'isotropicModes/whatI.dat'
  call read_1d_ascii(wI,trim(fname))
  
  write(fname,'(A)')trim(datadir)//'/'//'isotropicModes/kx.dat'
  call read_1d_ascii(k1,trim(fname))
  write(fname,'(A)')trim(datadir)//'/'//'isotropicModes/ky.dat'
  call read_1d_ascii(k2,trim(fname))
  write(fname,'(A)')trim(datadir)//'/'//'isotropicModes/kz.dat'
  call read_1d_ascii(k3,trim(fname))
  
  write(fname,'(A)')trim(datadir)//'/'//'isotropicModes/gmxloc.dat'
  call read_1d_ascii(xloc,trim(fname))
  write(fname,'(A)')trim(datadir)//'/'//'isotropicModes/gmyloc.dat'
  call read_1d_ascii(yloc,trim(fname))
  write(fname,'(A)')trim(datadir)//'/'//'isotropicModes/gmzloc.dat'
  call read_1d_ascii(zloc,trim(fname))

  ! Copy isotropic data to mode arrays
  call assert(size(k3) == 25600)

  uhatR = uR
  uhatI = uI
  vhatR = vR
  vhatI = vI
  whatR = wR
  whatI = wI

  kx = k1
  ky = k2
  kz = k3

  gmxloc = xloc
  gmyloc = yloc
  gmzloc = zloc

  deallocate(uR,uI,vR,vI,wR,wI,k1,k2,k3,xloc,yloc,zloc)

  ! Strain the modes
  call strainIsotropicModes()
  
  ! Read strained mode data
  write(fname,'(A)')trim(datadir)//'/'//'strainedModes/uhatR.dat'
  call read_1d_ascii(uR,trim(fname))
  write(fname,'(A)')trim(datadir)//'/'//'strainedModes/uhatI.dat'
  call read_1d_ascii(uI,trim(fname))
  write(fname,'(A)')trim(datadir)//'/'//'strainedModes/vhatR.dat'
  call read_1d_ascii(vR,trim(fname))
  write(fname,'(A)')trim(datadir)//'/'//'strainedModes/vhatI.dat'
  call read_1d_ascii(vI,trim(fname))
  write(fname,'(A)')trim(datadir)//'/'//'strainedModes/whatR.dat'
  call read_1d_ascii(wR,trim(fname))
  write(fname,'(A)')trim(datadir)//'/'//'strainedModes/whatI.dat'
  call read_1d_ascii(wI,trim(fname))
  
  write(fname,'(A)')trim(datadir)//'/'//'strainedModes/kx.dat'
  call read_1d_ascii(k1,trim(fname))
  write(fname,'(A)')trim(datadir)//'/'//'strainedModes/ky.dat'
  call read_1d_ascii(k2,trim(fname))
  write(fname,'(A)')trim(datadir)//'/'//'strainedModes/kz.dat'
  call read_1d_ascii(k3,trim(fname))
  
  write(fname,'(A)')trim(datadir)//'/'//'strainedModes/gmxloc.dat'
  call read_1d_ascii(xloc,trim(fname))
  write(fname,'(A)')trim(datadir)//'/'//'strainedModes/gmyloc.dat'
  call read_1d_ascii(yloc,trim(fname))
  write(fname,'(A)')trim(datadir)//'/'//'strainedModes/gmzloc.dat'
  call read_1d_ascii(zloc,trim(fname))
  
  ! Compare
  call assert(maxval(abs(uR - uhatR)) < tol, 'uhatR error')
  call assert(maxval(abs(uI - uhatI)) < tol, 'uhatI error')
  call assert(maxval(abs(vR - vhatR)) < tol, 'vhatR error')
  call assert(maxval(abs(vI - vhatI)) < tol, 'vhatI error')
  call assert(maxval(abs(wR - whatR)) < tol, 'whatR error')
  call assert(maxval(abs(wI - whatI)) < tol, 'whatI error')
  
  call assert(maxval(abs(k1 - kx)) < tol, 'whatI error')
  call assert(maxval(abs(k2 - ky)) < tol, 'whatI error')
  call assert(maxval(abs(k3 - kz)) < tol, 'whatI error')
  
  call assert(maxval(abs(xloc - gmxloc)) < tol, 'whatI error')
  call assert(maxval(abs(yloc - gmyloc)) < tol, 'whatI error')
  call assert(maxval(abs(zloc - gmzloc)) < tol, 'whatI error')
 
  call message('Test PASSED!') 
  call finalizeProblem()
end program
