include "../problems/gabor/turbHalfChannel_files/initialize.F90"

program test_gaborMode_largeScales
  use kind_parameters, only: rkind, clen
  use decomp_2d
  use mpi
  use domainSetup, only: nxLES, nyLES, nzLES, &
    gpLES, getStartAndEndIndices, dxLES, dyLES, dzLES, periodic
  use fortran_assert, only: assert
  use basic_io, only: read_1d_ascii
  use hdf5 
  use hdf5_fcns, only: read_h5_chunk_data
  use gaborIO_mod, only: readLargeScaleData
  use DerivativesMod, only: derivatives

  implicit none

  integer :: ierr, ioUnit
  real(rkind), dimension(:,:,:), allocatable :: Uascii, Vascii, Wascii
  real(rkind), dimension(:,:,:,:,:), allocatable :: gradUascii
  real(rkind), dimension(:), allocatable :: Uascii1, Vascii1, Wascii1
  real(rkind), dimension(:), allocatable :: gradUascii1
  real(rkind), dimension(:,:,:), allocatable :: U, V, W
  real(rkind), dimension(:,:,:,:,:), allocatable :: gradU
  character(len=clen) :: fname, datadir, inputfile
  integer :: ist, ien, jst, jen, kst, ken

  ! Derivatives
  type(derivatives) :: grad
  character(len=4) :: method_x = 'four'
  character(len=4) :: method_y = 'four'
  character(len=4) :: method_z = 'cd06'
  
  ! HDF5 variables
  integer(HID_T) :: memtype
  integer(HSIZE_T), dimension(3) :: dimsf, chunk_dims, stride, count1, block, offset
  character(len=2) :: dsetname
  integer :: dset_rank

  namelist /IO/ datadir

  ! Initialize MPI, HDF5, & generate numerical mesh
  call initializeProblem(inputfile)
  !call initializeExternalLibraries(ierr)
  !call MPI_Init(ierr)

  ! Initialize HDF5
  !call H5open_f(ierr)

  ! Get the input file path and file name
  !call GETARG(1,inputfile)
  
  ! Setup the domain
  !call setupDomainXYperiodic(inputfile)
  call getStartAndEndIndices(gpLES,ist,ien,jst,jen,kst,ken)

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
    
    ! Read partitioned fields
    allocate(U(ist:ien,jst:jen,kst:ken))
    allocate(V(ist:ien,jst:jen,kst:ken))
    allocate(W(ist:ien,jst:jen,kst:ken))
    
    write(fname,'(A)')trim(datadir)//'/'//'LargeScaleVelocity.h5'
    call readLargeScaleData(trim(fname),U,V,W,nxLES,nyLES,nzLES,gpLES)
    !U = Uascii(ist:ien,jst:jen,kst:ken) 
    !V = Vascii(ist:ien,jst:jen,kst:ken) 
    !W = Wascii(ist:ien,jst:jen,kst:ken)

  ! First, compute velocity gradient
  ! NOTE: gradU is stored as dUidXj = gradU(i,j)
  allocate(gradU(3,3,ist:ien,jst:jen,kst:ken))

  ! Initialize derivative object
  call grad%init(gpLES,dxLES,dyLES,dzLES,periodic(1),periodic(2),periodic(3),&
    method_x, method_y, method_z)

  ! Compute velocity gradient
  call grad%ddx(U,gradU(1,1,:,:,:))
  call grad%ddx(V,gradU(2,1,:,:,:))
  call grad%ddx(W,gradU(3,1,:,:,:))

  call grad%ddy(U,gradU(1,2,:,:,:))
  call grad%ddy(V,gradU(2,2,:,:,:))
  call grad%ddy(W,gradU(3,2,:,:,:))

  call grad%ddz(U,gradU(1,3,:,:,:))
  call grad%ddz(V,gradU(2,3,:,:,:))
  call grad%ddz(W,gradU(3,3,:,:,:))

  ! Verify this is the same as MATLAB
  print*, nrank, jst, kst
  print*, nrank, "Ux", maxval(abs(gradU(1,1,:,:,:) - gradUascii(1,1,ist:ien,jst:jen,kst:ken)))
  print*, nrank, "Vx", maxval(abs(gradU(2,1,:,:,:) - gradUascii(2,1,ist:ien,jst:jen,kst:ken)))
  print*, nrank, "Wx", maxval(abs(gradU(3,1,:,:,:) - gradUascii(3,1,ist:ien,jst:jen,kst:ken)))
  print*, nrank, "Uy", maxval(abs(gradU(1,2,:,:,:) - gradUascii(1,2,ist:ien,jst:jen,kst:ken)))
  print*, nrank, "Vy", maxval(abs(gradU(2,2,:,:,:) - gradUascii(2,2,ist:ien,jst:jen,kst:ken)))
  print*, nrank, "Wy", maxval(abs(gradU(3,2,:,:,:) - gradUascii(3,2,ist:ien,jst:jen,kst:ken)))
  print*, nrank, "Uz", maxval(abs(gradU(1,3,:,:,:) - gradUascii(1,3,ist:ien,jst:jen,kst:ken)))
  print*, nrank, "Vz", maxval(abs(gradU(2,3,:,:,:) - gradUascii(2,3,ist:ien,jst:jen,kst:ken)))
  print*, nrank, "Wz", maxval(abs(gradU(3,3,:,:,:) - gradUascii(3,3,ist:ien,jst:jen,kst:ken)))
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  call assert(maxval(abs(gradU - gradUascii(:,:,ist:ien,jst:jen,kst:ken))) < 1.d-11,'Test FAILED')
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  if (nrank == 0) print*, "Test PASSED!"
  call MPI_Barrier(MPI_COMM_WORLD,ierr)

  ! Free up memory
  deallocate(Uascii,Uascii1,Vascii,Vascii1,Wascii,Wascii1)
  deallocate(gradUascii1,gradUascii)
  deallocate(U,V,W,gradU) 
  !call finalizeDomainSetup()
  !call H5close_f(ierr)
  !call MPI_Finalize(ierr)
  call finalizeProblem()
end program
