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
  use gaborIO_mod, only: readVelocityField, writeFields
  use DerivativesMod, only: derivatives

  implicit none

  integer :: ierr, ioUnit
  real(rkind), dimension(:,:,:), allocatable :: Uascii, Vascii, Wascii
  real(rkind), dimension(:,:,:,:,:), allocatable :: gradUascii
  real(rkind), dimension(:), allocatable :: Uascii1, Vascii1, Wascii1
  real(rkind), dimension(:), allocatable :: gradUascii1
  real(rkind), dimension(:,:,:), allocatable :: U, V, W
  real(rkind), dimension(:,:,:), allocatable :: UinY
  real(rkind), dimension(:,:,:), allocatable :: UinZ
  real(rkind), dimension(:,:,:), allocatable :: dUinY, dUinZ
  real(rkind), dimension(:,:,:,:,:), allocatable :: gradU
  character(len=clen) :: fname, datadir, inputfile
  integer :: ist, ien, jst, jen, kst, ken, isz, jsz, ksz

  ! Derivatives
  type(derivatives) :: grad
  character(len=4) :: method_x = 'four'
  character(len=4) :: method_y = 'four'
  character(len=4) :: method_z = 'cd06'

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
    
    ! Read partitioned fields
    allocate(U(isz,jsz,ksz))
    allocate(V(isz,jsz,ksz))
    allocate(W(isz,jsz,ksz))
    
    write(fname,'(A)')trim(datadir)//'/'//'LargeScaleVelocity.h5'
    call readVelocityField(trim(fname),U,V,W,gpLES)

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

  ! First allocate memory and transpose arrays to use operators in specified
  ! direction
  allocate(UinY(gpLES%ysz(1),gpLES%ysz(2),gpLES%ysz(3)))
  allocate(dUinY(gpLES%ysz(1),gpLES%ysz(2),gpLES%ysz(3)))
  call transpose_x_to_y(U,UinY,gpLES)
  call grad%ddy(UinY,dUinY)
  call transpose_y_to_x(dUinY,gradU(1,2,:,:,:),gpLES)
  
  call transpose_x_to_y(V,UinY,gpLES)
  call grad%ddy(UinY,dUinY)
  call transpose_y_to_x(dUinY,gradU(2,2,:,:,:),gpLES)
  
  call transpose_x_to_y(W,UinY,gpLES)
  call grad%ddy(UinY,dUinY)
  call transpose_y_to_x(dUinY,gradU(3,2,:,:,:),gpLES)

  allocate(UinZ(gpLES%zsz(1),gpLES%zsz(2),gpLES%zsz(3)))
  allocate(dUinZ(gpLES%zsz(1),gpLES%zsz(2),gpLES%zsz(3)))
  call transpose_x_to_y(U,UinY,gpLES)
  call transpose_y_to_z(UinY,UinZ,gpLES)
  call grad%ddz(UinZ,dUinZ)
  call transpose_z_to_y(dUinZ,dUinY,gpLES)
  call transpose_y_to_x(dUinY,gradU(1,3,:,:,:),gpLES)

  call transpose_x_to_y(V,UinY,gpLES)
  call transpose_y_to_z(UinY,UinZ,gpLES)
  call grad%ddz(UinZ,dUinZ)
  call transpose_z_to_y(dUinZ,dUinY,gpLES)
  call transpose_y_to_x(dUinY,gradU(2,3,:,:,:),gpLES)

  call transpose_x_to_y(W,UinY,gpLES)
  call transpose_y_to_z(UinY,UinZ,gpLES)
  call grad%ddz(UinZ,dUinZ)
  call transpose_z_to_y(dUinZ,dUinY,gpLES)
  call transpose_y_to_x(dUinY,gradU(3,3,:,:,:),gpLES)

  ! Verify this is the same as MATLAB
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  call assert(maxval(abs(gradU - gradUascii(:,:,ist:ien,jst:jen,kst:ken))) < 1.d-11,'Test FAILED')
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  if (nrank == 0) print*, "Test PASSED!"
  call MPI_Barrier(MPI_COMM_WORLD,ierr)

  ! Write gradient to disk
  write(fname,'(A)')trim(datadir)//'/'//'LargeScaleGradientOUT.h5'
  call writeFields(trim(fname),gradU,'/gradU',gpLES)

  ! Free up memory
  deallocate(Uascii,Uascii1,Vascii,Vascii1,Wascii,Wascii1)
  deallocate(gradUascii1,gradUascii)
  deallocate(U,V,W,gradU)
  deallocate(UinY,dUinY) 
  deallocate(UinZ,dUinZ) 
  call finalizeProblem()
end program
