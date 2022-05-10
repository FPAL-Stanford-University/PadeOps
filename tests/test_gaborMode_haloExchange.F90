program test_gaborMode_haloExchange
  use kind_parameters, only: rkind, clen
  use decomp_2d
  use decomp_2d_io
  use mpi
  use GaborModeRoutines, only: haloExchange
  use domainSetup, only: getStartAndEndIndices, periodic, decomp2Dpencil,&
    nysupp, nxsupp, nzsupp
  use largeScalesMod, only: writeFields
  use fortran_assert, only: assert
  use hdf5
  use exits, only: message
  implicit none

  real(rkind), dimension(:,:,:), allocatable :: A, Aout
  real(rkind), dimension(:,:,:), allocatable :: Ahalo
  integer :: ierr
  type(decomp_info) :: gp, gph
  integer :: ist, ien, jst, jen, kst, ken, isz, jsz, ksz
  character(len=clen) :: fname1, fname2, fname3, fname4, fname5
  character(len=clen) :: outputdir
  integer :: nx, ny, nz, nprocX, nprocY, nprocZ

  call MPI_Init(ierr)
  call H5open_f(ierr)

  ! IO stuff
  outputdir = '/work2/06632/ryanhass/stampede2/Enrichment/GaborKS_V2Data/tests/MPI/'
  fname1 = trim(outputdir)//'testMPIhaloExchangeAout.h5'
  fname2 = trim(outputdir)//'testMPIhaloExchangeAout.out'
  fname3 = trim(outputdir)//'testMPIhaloExchangeAin.h5'
  fname4 = trim(outputdir)//'testMPIhaloExchangeAin.out'
  fname5 = trim(outputdir)//'testMPIhaloExchangeAinBeforeExchange.h5'
  
  periodic = [.true.,.true.,.false.]
  nx = 64
  ny = 64
  nz = 32
  call decomp_2d_init(nx,ny,nz,2,2,periodic)
  call get_decomp_info(gp)
  call getStartAndEndIndices(gp,ist,ien,jst,jen,kst,ken,isz,jsz,ksz)
  allocate(Aout(ist,jsz,ksz))

  ! Domain setup
  decomp2Dpencil = 'x';
  nxsupp = 8;
  nysupp = 8;
  nzsupp = 8;

 !TODO: need new grid partition for Ah to match uFh definitions
  nprocX = 1;
  nprocY = 2;
  nprocZ = 2;
  call decomp_info_init(nx + nprocX*nxsupp, ny+nprocY*nysupp, nz+(nprocZ-1)*nzsupp,gph)
  call getStartAndEndIndices(gph,ist,ien,jst,jen,kst,ken,isz,jsz,ksz)
  allocate(A(isz,jsz,ksz))
  A = nrank + 1

  call assert(ist+nxsupp < ien,'ist+nxsupp < ien')
  call assert(jst+nysupp < jen,'jst+nysupp < jen')
  call assert(kst+nzsupp < ken,'kst+nzsupp < ken')

  call update_halo(A,Ahalo,level=nxsupp/2,opt_decomp=gph,opt_global=.true.)
  call writeFields(trim(fname5),A,'Ai',gph)
  
  call haloExchange(A,Ahalo)
 
  if (kst == 1) then
    Aout = A(1+nxsupp/2:isz-nxsupp/2,1+nysupp/2:jsz-nysupp/2,&
      1:ksz-nzsupp/2)
  else
    Aout = A(1+nxsupp/2:isz-nxsupp/2,1+nysupp/2:jsz-nysupp/2,&
      1+nzsupp/2:ksz)
  endif
  
  call writeFields(trim(fname1),Aout,'Ao',gp)
  call decomp_2d_write_one(1,Aout,trim(fname2),gp)

  call writeFields(trim(fname3),A,'Ai',gph)
  call decomp_2d_write_one(1,A,trim(fname4),gph)
  deallocate(A, Ahalo)
  call decomp_info_finalize(gp)
  call decomp_2d_finalize
  call H5close_f(ierr)
  call MPI_Finalize(ierr)
end program
