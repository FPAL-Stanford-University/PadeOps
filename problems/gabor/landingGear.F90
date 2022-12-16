#include "landingGear_files/initialize.F90"

program landingGear
  use kind_parameters,         only: clen, rkind
  use enrichmentMod,           only: enrichmentOperator, nthreads, interpAndAddGrids 
  use IncompressibleGrid,      only: igrid
  use auxiliary_openmp_subs,   only: GetArguments
  use mpi
  use interpolatorMod,         only: interpolator
  use reductions, only: p_sum
  implicit none

  character(len=clen) :: inputS1, inputS2, inputS3, inputG12, inputG23
  integer :: ierr, provided
  type(igrid) :: S1, S2, S3, S4
  type(enrichmentOperator) :: enrich12, enrich23
  type(interpolator) :: interpS1ToS2, interpS2ToS3
  real(rkind), dimension(:,:,:), allocatable :: u, v, w
  integer :: tid = 0  
  logical :: dumpIndividualScales = .true. 
  call MPI_Init_thread(MPI_THREAD_FUNNELED,provided,ierr)
  
  call GETARG(1,inputS1)
  call GETARG(2,inputS2)
  call GETARG(3,inputS3)
  call GETARG(4,inputG12)
  call GETARG(5,inputG23)
  call GetArguments(nthreads,6)

  call S1%init(inputS1, .true.) 
  
  call S2%init(inputS2, .false.)
  !call S3%init(inputS3, .false.)
    
  
  call interpS1ToS2%init(S1%gpC, S2%gpC, S1%x1d, S1%y1d, S1%z1d, & 
        S2%x1d, S2%y1d, S2%z1d)
  !call interpS2ToS3%init(S2%gpC, S3%gpC, S2%x1d, S2%y1d, S2%z1d, & 
  !      S3%x1d, S3%y1d, S3%z1d)

  call enrich12%init(S2,S1,inputG12)
  !call enrich23%init(S3,S2,inputG23)
 
  ! Scale 1 : read from disk 
  call S1%initLargeScales(tid,readGradients=.true.) 
  !call S1%ComputeCD06Gradients([.false.,.false.,.false.])
  
  ! Scale 2 : Generate using modes 
  call enrich12%generateModes()
  call enrich12%renderVelocity()
  call S2%ComputeCD06Gradients([.false.,.false.,.false.])
  if (dumpIndividualScales) call enrich12%dumpSmallScales()

  print*, p_sum(sum(abs(S2%u)))
  !! Update Scale 2 for Scale 3
  !call interpAndAddGrids(S1, S2, interpS1ToS2)

  !! Scale 3: Generate Using modes 
  !call enrich23%generateModes()
  !call enrich23%renderVelocity()
  !call S3%ComputeCD06Gradients([.false.,.false.,.false.])
  !if (dumpIndividualScales) call enrich23%dumpSmallScales()


  call enrich12%destroy()
  !call enrich23%destroy()
  call S1%destroy()
  call S2%destroy()
  
  !call S3%destroy()


  call MPI_Finalize(ierr)

end program landingGear
