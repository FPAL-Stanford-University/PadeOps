program test_gaborMode_haloExchangeCustom
  use mpi
  use decomp_2d
  use kind_parameters, only: rkind, clen
  use domainSetup, only: decomp2Dpencil, getStartAndEndIndices, nxsupp, &
    nysupp, nzsupp, periodic
  use GaborModeRoutines, only: buff1, buff2, buff3, buff4, haloExchangeMPI, &
    istAll, ienAll, jstAll, jenAll, kstAll, kenAll, finalizeGaborModes, &
    sendReqLo, sendReqHi, recvReqLo, recvReqHi
  use largeScalesMod, only: writeFields
  use hdf5
  use fortran_assert, only: assert
  implicit none

  real(rkind), dimension(:,:,:), allocatable :: uh, u
  integer :: nx, ny, nz
  type(decomp_info) :: gp
  integer :: ierr, ist, ien, jst, jen, kst, ken, isz, jsz, ksz
  integer :: ksth, kenh, kszh
  integer, dimension(:), allocatable :: myIst, myIen, myJst, myJen, myKst, myKen
  character(len=clen) :: fname

  call MPI_Init(ierr)
  call H5open_f(ierr)

  ! Output file name
  fname = '/work2/06632/ryanhass/stampede2/Enrichment/GaborKS_V2Data/tests/MPI/testMPIhaloExchange_4PE_YZ.h5'
  
  ! Domain setup
  nx = 32
  ny = 32
  nz = 64

  nxsupp = 4
  nysupp = 4
  nzsupp = 4

  periodic = [.true.,.true.,.false.]

  ! Initialize domain decomposition for MPI
  call decomp_2d_init(nx,ny,nz+1,2,2,periodic)
  call get_decomp_info(gp)

  ! X-decomposition
  decomp2Dpencil = 'x'
  call getStartAndEndIndices(gp,ist,ien,jst,jen,kst,ken,isz,jsz,ksz)
  
  ! Bounds check to make sure we don't have too many MPI processes for the data
  ! size
  call assert(ist+nxsupp < ien,'ist+nxsupp < ien -- domainSetup.F90')
  call assert(jst+nysupp < jen,'jst+nysupp < jen -- domainSetup.F90')
  call assert(kst+nzsupp < ken,'kst+nzsupp < ken -- domainSetup.F90')
  
  ! Allocate MPI communication arrays
  allocate(istAll(nproc),ienAll(nproc))
  allocate(jstAll(nproc),jenAll(nproc))
  allocate(kstAll(nproc),kenAll(nproc))

  allocate(myIst(nproc),myIen(nproc))
  allocate(myJst(nproc),myJen(nproc))
  allocate(myKst(nproc),myKen(nproc))
  
  myIst = ist; myIen = ien
  myJst = jst; myJen = jen
  myKst = kst; myKen = ken

  ! Send the start and end indices to all other processes
  call MPI_Alltoall(myIst,1,MPI_INTEGER,istAll,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
  call MPI_Alltoall(myIen,1,MPI_INTEGER,ienAll,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
  call MPI_Alltoall(myJst,1,MPI_INTEGER,jstAll,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
  call MPI_Alltoall(myJen,1,MPI_INTEGER,jenAll,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
  call MPI_Alltoall(myKst,1,MPI_INTEGER,kstAll,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
  call MPI_Alltoall(myKen,1,MPI_INTEGER,kenAll,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)

  deallocate(myIst,myIen,myJst,myJen,myKst,myKen)

  ! Output array
  allocate(u(isz,jsz,ksz))

  ! Halo-padded array
  ksth = max(1,kst-nzsupp/2)
  kenh = min(nz+1,ken+nzsupp/2)
  kszh = kenh - ksth + 1
  allocate(uh(isz+nxsupp,jsz+nysupp,kszh))
  allocate(buff1(isz+nxsupp,jsz+nysupp,nzsupp/2))
  allocate(buff2(isz+nxsupp,jsz+nysupp,nzsupp/2))
  allocate(buff3(isz+nxsupp,nysupp/2,kszh))
  allocate(buff4(isz+nxsupp,nysupp/2,kszh))

  allocate(sendReqLo(kszh), sendReqHi(kszh))
  allocate(recvReqLo(kszh), recvReqHi(kszh))

  uh = real(nrank+1,rkind)
  call haloExchangeMPI(uh)

  if (kszh == ksz) then
    u = uh(1+nxsupp/2:isz+nxsupp/2,1+nysupp/2:jsz+nysupp/2,1:ksz)
  elseif (ksth == 1) then
    u = uh(1+nxsupp/2:isz+nxsupp/2,1+nysupp/2:jsz+nysupp/2,1:ksz)
  else
    u = uh(1+nxsupp/2:isz+nxsupp/2,1+nysupp/2:jsz+nysupp/2,1+nzsupp/2:ksz+nzsupp/2)
  end if
  
  call writeFields(trim(fname),u,'/u',gp)
  deallocate(uh,u)
  call finalizeGaborModes()
  call decomp_info_finalize(gp)
  call H5close_f(ierr)
  call MPI_Finalize(ierr)
end program
