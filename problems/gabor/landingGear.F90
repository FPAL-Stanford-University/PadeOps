#include "landingGear_files/initialize.F90"

program landingGear
  use kind_parameters,         only: clen, rkind
  use enrichmentMod,           only: enrichmentOperator, nthreads
  use IncompressibleGrid,      only: igrid
  use auxiliary_openmp_subs,   only: GetArguments
  use mpi
  use interpolatorMod,         only: interpolator

  implicit none

  character(len=clen) :: inputS1, inputS2, inputS3, inputG12, inputG23
  integer :: ierr, provided
  type(igrid) :: S1, S2, S3
  type(enrichmentOperator) :: enrich12, enrich23
  type(interpolator) :: interp
  real(rkind), dimension(:,:,:), allocatable :: u, v, w
      
  call MPI_Init_thread(MPI_THREAD_FUNNELED,provided,ierr)
  
  call GETARG(1,inputS1)
  call GETARG(2,inputS2)
  call GETARG(3,inputG12)
  call GetArguments(nthreads,4)

  call S1%init(inputS1, .true.) 
  call S2%init(inputS2, .false.)
  call S3%init(inputS3, .false.)
    
  
  call enrich12%init(S2,S1,inputG12)
  call enrich12%renderVelocity()
  call enrich12%dumpSmallScales()
 
  ! Interpolate S1 to S2 grid and add the fields 
  allocate(u(S2%gpC%xsz(1), S2%gpC%xsz(2), S2%gpC%xsz(3)))
  allocate(v(S2%gpC%xsz(1), S2%gpC%xsz(2), S2%gpC%xsz(3)))
  allocate(w(S2%gpC%xsz(1), S2%gpC%xsz(2), S2%gpC%xsz(3)))

  call interp%init(S1%gpC, S2%gpC, S1%mesh(:,1,1,1), S1%mesh(1,:,1,2), &
    S1%mesh(1,1,:,3), S2%mesh(:,1,1,1), S2%mesh(1,:,1,2), S2%mesh(1,1,:,3)) 
  call interp%LinInterp3D(S1%u ,u)
  call interp%LinInterp3D(S1%v ,v)
  call interp%LinInterp3D(S1%wC,w)
  S2%u  = S2%u  + u
  S2%v  = S2%v  + v
  S2%wC = S2%wC + w

  ! Do level three enrichment
  call enrich23%init(S3,S2,inputG23)
  call enrich23%renderVelocity()
  call enrich23%dumpSmallScales()

  call enrich12%destroy()
  call enrich23%destroy()
  call S1%destroy()
  call S2%destroy()
  call S3%destroy()

  deallocate(u,v,w)

  call MPI_Finalize(ierr)
end program landingGear
