#include "../problems/gabor/hitEnrich_files/initialize.F90"
program test_gaborMode_renderVelMPI
  use kind_parameters,       only: clen, rkind
  use enrichmentMod,         only: enrichmentOperator, nthreads, xDom, yDom, zDom
  use IncompressibleGrid,    only: igrid
  use auxiliary_openmp_subs, only: GetArguments
  use mpi
  use exits,                 only: message
  use decomp_2d,             only: nrank
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

  enrich%uhatR = 1.d0
  enrich%uhatI = 1.d0
  enrich%vhatR = 1.d0
  enrich%vhatI = 1.d0
  enrich%whatR = 1.d0
  enrich%whatI = 1.d0

  enrich%kx = 1.d0
  enrich%ky = 1.d0
  enrich%kz = 1.d0
  
  Lx = xDom(2) - xDom(1)
  Ly = yDom(2) - yDom(1)
  Lz = zDom(2) - zDom(1)
  
  enrich%x = 0.1d0
  enrich%y = 0.1d0
  enrich%z = 0.1d0

  !enrich%smallScales%u  = 0.d0
  !enrich%smallScales%v  = 0.d0
  !enrich%smallScales%wC = 0.d0
  !enrich%utmp = 0.d0
  !enrich%vtmp = 0.d0
  !enrich%wtmp = 0.d0
  !call enrich%renderLocalVelocity(enrich%x,enrich%y,enrich%z,enrich%kx,enrich%ky,enrich%kz,&
  !    enrich%uhatR,enrich%uhatI,enrich%vhatR,enrich%vhatI,enrich%whatR,enrich%whatI)
  call enrich%renderVelocity()
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  print*, nrank, "here"
  call enrich%dumpSmallScales()

  call enrich%destroy()
  call largeScales%destroy()
  call smallScales%destroy()

  call MPI_Finalize(ierr)
end program 
