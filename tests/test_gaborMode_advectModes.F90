#include "../problems/gabor/hitEnrich_files/initialize.F90"
program hitEnrich
  use kind_parameters,         only: clen, rkind
  use enrichmentMod,           only: enrichmentOperator, nthreads, zDom 
  use IncompressibleGrid,      only: igrid
  use auxiliary_openmp_subs,   only: GetArguments
  use mpi
  use reductions
  use decomp_2d,               only: nrank, nproc
  use constants,               only: two, pi
  use exits,                   only: message
  use fortran_assert,          only: assert
  implicit none

  character(len=clen) :: inputfileLS, inputfileSS, inputfileGM 
  integer :: ierr, provided
  type(igrid) :: largeScales, smallScales
  type(enrichmentOperator) :: enrich
  real(rkind), dimension(:,:,:), allocatable :: u1, v1, w1
  real(rkind) :: tol = 1.d-13
  real(rkind) :: Lz
      
  call MPI_Init_thread(MPI_THREAD_FUNNELED,provided,ierr)
  
  call GETARG(1,inputfileLS)
  call GETARG(2,inputfileSS)
  call GETARG(3,inputfileGM)
  call GetArguments(nthreads,4)

  call largeScales%init(inputfileLS, .true.) 
  call smallScales%init(inputfileSS, .false.)
  call enrich%init(smallScales,largeScales,inputfileGM)
  Lz = zDom(2) - zDom(1)

  ! Set uniform flow in one coordinate direction
  largeScales%v       = two*pi !enrich%smallScales%dy/enrich%dt 
  largeScales%u       = 0.d0
  largeScales%wC      = 0.d0
  largeScales%duidxjC = 0.d0
  
  enrich%tid = 0
  call enrich%wrapupTimeStep()
  
  allocate(u1(size(enrich%smallScales%u, 1),size(enrich%smallScales%u, 2),size(enrich%smallScales%u, 3)))
  allocate(v1(size(enrich%smallScales%v, 1),size(enrich%smallScales%v, 2),size(enrich%smallScales%v, 3)))
  allocate(w1(size(enrich%smallScales%wC,1),size(enrich%smallScales%wC,2),size(enrich%smallScales%wC,3)))
  u1 = enrich%smallScales%u
  v1 = enrich%smallScales%v
  w1 = enrich%smallScales%wC
  
  do while (enrich%continueSimulation())
    call enrich%updateLargeScales(timeAdvance=.false.,initializing=.false.)
    call enrich%advanceTime()
    call enrich%wrapupTimeStep()

    call assert(enrich%nmodes == size(enrich%uhatR),'size mismatch')

    if (mod(enrich%tid,10) == 0) then
      call message('step ',enrich%tid,' of ',enrich%tidStop)
    end if
  end do

  print*, "max u difference:", p_maxval(maxval(abs(enrich%smallScales%u - u1)))
  print*, "max v difference:", p_maxval(maxval(abs(enrich%smallScales%v - v1)))
  print*, "max w difference:", p_maxval(maxval(abs(enrich%smallScales%wC - w1)))
  call assert(p_maxval(maxval(abs(enrich%smallScales%u - u1))) < tol)
  call assert(p_maxval(maxval(abs(enrich%smallScales%v - v1))) < tol)
  call assert(p_maxval(maxval(abs(enrich%smallScales%wC - w1))) < tol)
  call message("Test PASSED!")
  deallocate(u1,v1,w1)
  call enrich%destroy()
  call largeScales%destroy()
  call smallScales%destroy()

  call MPI_Finalize(ierr)
end program hitEnrich
