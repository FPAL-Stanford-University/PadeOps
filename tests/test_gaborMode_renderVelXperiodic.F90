#include "../problems/gabor/hitEnrich_files/initialize.F90"
program test_gaborMode_renderVelXperiodic
  use kind_parameters,         only: clen, rkind
  use enrichmentMod,           only: enrichmentOperator, nthreads, xDom, yDom, zDom
  use IncompressibleGrid,      only: igrid
  use auxiliary_openmp_subs,   only: GetArguments
  use mpi
  use exits,                   only: message
  implicit none

  character(len=clen) :: inputfileLS, inputfileSS, inputfileGM 
  character(len=clen) :: datadir, fname, outputdir
  integer :: ioUnit, ierr, provided
  type(igrid) :: largeScales, smallScales
  type(enrichmentOperator) :: enrich
  real(rkind) :: Lx, Ly, Lz, xWindow
      
  call MPI_Init_thread(MPI_THREAD_FUNNELED,provided,ierr)
  
  call GETARG(1,inputfileLS)
  call GETARG(2,inputfileSS)
  call GETARG(3,inputfileGM)
  call GetArguments(nthreads)
  
  call largeScales%init(inputfileLS, .true.) 
  call smallScales%init(inputfileSS, .false.)

  call enrich%init(smallScales,largeScales,inputfileGM)

  enrich%uhatR = 0.d0
  enrich%uhatI = 0.d0
  enrich%vhatR = 0.d0
  enrich%vhatI = 0.d0
  enrich%whatR = 0.d0
  enrich%whatI = 0.d0

  enrich%uhatR(1) = 1.d0
  enrich%uhatI(1) = 1.d0
  enrich%vhatR(1) = 1.d0
  enrich%vhatI(1) = 1.d0
  enrich%whatR(1) = 1.d0
  enrich%whatI(1) = 1.d0

  enrich%kx(1) = 1.d0
  enrich%ky(1) = 0.d0
  enrich%kz(1) = 0.d0

  LX = xDom(2) - xDom(1)
  Ly = yDom(2) - yDom(1)
  Lz = zDom(2) - zDom(1)
  xWindow = enrich%nxsupp*smallScales%dx

  enrich%x(1) = Lx-xWindow/8.d0!xWindow/8.d0
  enrich%y(1) = Ly/2.d0
  enrich%z(1) = Lz/2.d0
    
  call enrich%renderVelocity()
  call enrich%dumpSmallScales()

  call enrich%destroy()
  call largeScales%destroy()
  call smallScales%destroy()

  call MPI_Finalize(ierr)
end program 
