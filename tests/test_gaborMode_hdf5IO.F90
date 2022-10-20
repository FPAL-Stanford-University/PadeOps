program test_gaborMode_hdf5IO
  use kind_parameters, only: rkind, clen
  use hdf5
  use hdf5_fcns, only: read_h5_chunk_data
  use mpi
  use domainSetup, only: setupDomainXYperiodic, finalizeDomainSetup, nxLES, nyLES, nzLES, &
    gpLES, getStartAndEndIndices
  use decomp_2d
  use basic_io, only: read_1d_ascii
  use fortran_assert, only: assert
  use gaborIO_mod, only: readFields, writeFields
  
  implicit none
  
  character(len=clen) :: fname, datadir, inputfile
  integer :: ierr, ioUnit
  integer :: ist, ien, jst, jen, kst, ken, isz, jsz, ksz
  real(rkind), dimension(:,:,:), allocatable :: Uascii, Vascii, Wascii
  real(rkind), dimension(:,:,:,:,:), allocatable :: gradUascii
  real(rkind), dimension(:), allocatable :: Uascii1, Vascii1, Wascii1, gradUascii1
  real(rkind), dimension(:,:,:), allocatable :: U, V, W
  real(rkind), dimension(:,:,:,:,:), allocatable :: gradU
  
  namelist /IO/ datadir

  ! Initialize MPI
  call MPI_Init(ierr)
  call H5open_f(ierr)

  ! Get the input file path and file name
  call GETARG(1,inputfile)
  
  ! Setup the domain
  call setupDomainXYperiodic(inputfile)

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

    allocate(gradUascii(3,3,nxLES,nyLES,nzLES))
    write(fname,'(A9)')'gradU.dat'
    call read_1d_ascii(gradUascii1,trim(datadir)//'/'//trim(fname))
    gradUascii = reshape(gradUascii1,[3,3,nxLES,nyLES,nzLES])
    
    ! Read partitioned fields
    call getStartAndEndIndices(gpLES,ist,ien,jst,jen,kst,ken,isz,jsz,ksz)
    allocate(U(isz,jsz,ksz))
    allocate(V(isz,jsz,ksz))
    allocate(W(isz,jsz,ksz))
    allocate(gradU(3,3,isz,jsz,ksz))

    ! HDF5 specific variables
    write(fname,'(A)') trim(datadir)//'/'//'LargeScaleVelocity.h5'
    call readFields(trim(fname),U,V,W,'/u','/v','/w',gpLES)

    ! Read in data
    call assert(maxval(abs(U-Uascii(ist:ien,jst:jen,kst:ken)))<1.d-14,'U diff')
    call assert(maxval(abs(V-Vascii(ist:ien,jst:jen,kst:ken)))<1.d-14,'V diff')
    call assert(maxval(abs(W-Wascii(ist:ien,jst:jen,kst:ken)))<1.d-14,'W diff')
    
    write(fname,'(A)') trim(datadir)//'/'//'LargeScaleGradient.h5'
    call readFields(trim(fname),gradU,'/gradU',gpLES)

    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    if (nrank == 0) print*, "Reading test PASSED!"
    call MPI_Barrier(MPI_COMM_WORLD,ierr)

  ! Run writing test
    write(fname,'(A)') trim(datadir)//'/'//'LargeScaleVelocityOUT.h5'
    call writeFields(trim(fname),U,V,W,'/u','/v','/w',gpLES)  

    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    if (nrank == 0) print*, "Finished writing data."
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
  ! Deallocate all memory
  deallocate(U,V,W,gradU)
  deallocate(Uascii1,Vascii1,Wascii1,gradUascii1)
  deallocate(Uascii,Vascii,Wascii,gradUascii)
  call finalizeDomainSetup()
  call H5close_f(ierr)
  call MPI_Finalize(ierr)
end program
