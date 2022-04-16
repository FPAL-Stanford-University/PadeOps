include "../problems/gabor/turbHalfChannel_files/initialize.F90"

program test_gaborMode_largeScales
  use kind_parameters, only: rkind, clen
  use decomp_2d
  use mpi
  use domainSetup, only: nxLES, nyLES, nzLES, gpLESb, gpLES, getStartAndEndIndices, &
    finalizeDomainSetup, setupDomainXYperiodic
  use fortran_assert, only: assert
  use basic_io, only: read_1d_ascii
  use hdf5 
  use hdf5_fcns, only: read_h5_chunk_data
  use gaborIO_mod, only: writeFields
  use DerivativesMod, only: derivatives
  use largeScalesMod, only: gradU, getLargeScaleData, computeLargeScaleGradient, &
    finalizeLargeScales, initLargeScales, U, V, W

  implicit none

  integer :: ierr, ioUnit
  real(rkind), dimension(:,:,:), allocatable :: Uascii, Vascii, Wascii
  real(rkind), dimension(:,:,:,:,:), allocatable :: gradUascii
  real(rkind), dimension(:), allocatable :: Uascii1, Vascii1, Wascii1
  real(rkind), dimension(:), allocatable :: gradUascii1
  character(len=clen) :: fname, datadir, inputfile
  integer :: ist, ien, jst, jen, kst, ken, isz, jsz, ksz

  namelist /IO/ datadir

  ! Initialize MPI, HDF5, & generate numerical mesh
  call initializeProblem(inputfile)
  
  call getStartAndEndIndices(gpLES,ist,ien,jst,jen,kst,ken,isz,jsz,ksz)

  ! Get ground truth from MATLAB
    ! Read inputfile
      open(unit=ioUnit, file=trim(inputfile), form='FORMATTED', iostat=ierr)
      read(unit=ioUnit, NML=IO)
      close(ioUnit)

    ! Read full fields
    allocate(Uascii(nxLES,nyLES,nzLES))
    allocate(Vascii(nxLES,nyLES,nzLES))
    allocate(Wascii(nxLES,nyLES,nzLES))
    
    write(fname,'(A5)')'U.dat'
    call read_1d_ascii(Uascii1,trim(datadir)//'/'//trim(fname))
    Uascii = reshape(Uascii1,[nxLES,nyLES,nzLES])
    write(fname,'(A5)')'V.dat'
    call read_1d_ascii(Vascii1,trim(datadir)//'/'//trim(fname))
    Vascii = reshape(Vascii1,[nxLES,nyLES,nzLES])
    write(fname,'(A5)')'W.dat'
    call read_1d_ascii(Wascii1,trim(datadir)//'/'//trim(fname))
    Wascii = reshape(Wascii1,[nxLES,nyLES,nzLES])


    !allocate(gradUascii(3,3,nxLES+1,nyLES+1,nzLES+2))
    allocate(gradUascii(3,3,nxLES,nyLES,nzLES))
    write(fname,'(A9)')'gradU.dat'
    call read_1d_ascii(gradUascii1,trim(datadir)//'/'//trim(fname))
    !gradUascii = reshape(gradUascii1,[3,3,nxLES+1,nyLES+1,nzLES+2])
    gradUascii = reshape(gradUascii1,[3,3,nxLES,nyLES,nzLES])
    
    write(fname,'(A)')trim(datadir)//'/'//'LargeScaleVelocity.h5'
    call getLargeScaleData(fname)

  ! First, compute velocity gradient
  call computeLargeScaleGradient()

  ! Verify this is the same as MATLAB
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  call assert(maxval(abs(gradU - gradUascii(:,:,ist:ien,jst:jen,kst:ken))) < 1.d-11,'Test FAILED')
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  if (nrank == 0) print*, "PART1 of test PASSED!"
  call MPI_Barrier(MPI_COMM_WORLD,ierr)

  ! Write gradient to disk
  write(fname,'(A)')trim(datadir)//'/'//'LargeScaleGradientOUT.h5'
  call writeFields(trim(fname),gradU,'/gradU',gpLES)

  ! Free up memory
  deallocate(Uascii,Uascii1,Vascii,Vascii1,Wascii,Wascii1)
  deallocate(gradUascii1,gradUascii)
  call finalizeDomainSetup()
  call finalizeLargeScales()
  !call finalizeProblem()

!! PART 2 OF TEST 
  ! Now make sure we can alternativel read the velocity and gradient from disk
  !call initializeProblem(inputfile)
  call setupDomainXYperiodic(inputfile)
  call initLargeScales(.false.)
  
  call getStartAndEndIndices(gpLESb,ist,ien,jst,jen,kst,ken,isz,jsz,ksz)

  ! Get ground truth from MATLAB
    ! Read inputfile
      open(unit=ioUnit, file=trim(inputfile), form='FORMATTED', iostat=ierr)
      read(unit=ioUnit, NML=IO)
      close(ioUnit)

    ! Read full fields
    allocate(Uascii(nxLES+1,nyLES+1,nzLES+2))
    allocate(Vascii(nxLES+1,nyLES+1,nzLES+2))
    allocate(Wascii(nxLES+1,nyLES+1,nzLES+2))
    
    write(fname,'(A6)')'UB.dat'
    call read_1d_ascii(Uascii1,trim(datadir)//'/'//trim(fname))
    Uascii = reshape(Uascii1,[nxLES+1,nyLES+1,nzLES+2])
    write(fname,'(A6)')'VB.dat'
    call read_1d_ascii(Vascii1,trim(datadir)//'/'//trim(fname))
    Vascii = reshape(Vascii1,[nxLES+1,nyLES+1,nzLES+2])
    write(fname,'(A6)')'WB.dat'
    call read_1d_ascii(Wascii1,trim(datadir)//'/'//trim(fname))
    Wascii = reshape(Wascii1,[nxLES+1,nyLES+1,nzLES+2])

    allocate(gradUascii(3,3,nxLES+1,nyLES+1,nzLES+2))
    write(fname,'(A10)')'gradUb.dat'
    call read_1d_ascii(gradUascii1,trim(datadir)//'/'//trim(fname))
    gradUascii = reshape(gradUascii1,[3,3,nxLES+1,nyLES+1,nzLES+2])
    
    write(fname,'(A)')trim(datadir)//'/'//'LargeScaleVelocityB.h5'
    call getLargeScaleData(fname,.true.)

  ! First, compute velocity gradient
    write(fname,'(A)')trim(datadir)//'/'//'LargeScaleGradientB.h5'
  call computeLargeScaleGradient(trim(fname))
  
  ! Verify this is the same as MATLAB
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  call assert(maxval(abs(gradU - gradUascii(:,:,ist:ien,jst:jen,kst:ken))) < 1.d-11,'Test FAILED: gradU discrepency')
  call assert(maxval(abs(U - Uascii(ist:ien,jst:jen,kst:ken))) < 1.d-11,'Test FAILED: U discrepency')
  call assert(maxval(abs(V - Vascii(ist:ien,jst:jen,kst:ken))) < 1.d-11,'Test FAILED: V discrepency')
  call assert(maxval(abs(W - Wascii(ist:ien,jst:jen,kst:ken))) < 1.d-11,'Test FAILED: W discrepency')
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  if (nrank == 0) print*, "PART2 of test PASSED!"
  call MPI_Barrier(MPI_COMM_WORLD,ierr)

  ! Free up memory
  deallocate(Uascii,Uascii1,Vascii,Vascii1,Wascii,Wascii1)
  deallocate(gradUascii1,gradUascii)
  call finalizeProblem()
end program
