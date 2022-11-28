#include "../problems/gabor/hitEnrich_files/initialize.F90"
program test_gaborMode_MPImodeExchange
  use kind_parameters,         only: clen, rkind
  use enrichmentMod,           only: enrichmentOperator, nthreads
  use IncompressibleGrid,      only: igrid
  use auxiliary_openmp_subs,   only: GetArguments
  use mpi
  use decomp_2d,               only: nrank, nproc
  use fortran_assert,          only: assert
  implicit none

  character(len=clen) :: inputfileLS, inputfileSS, inputfileGM 
  integer ::  ierr, provided
  type(igrid) :: largeScales, smallScales
  type(enrichmentOperator), target :: enrich
  real(rkind), dimension(:), allocatable :: x, y, z, kx, ky, kz, uR, uI, vR, vI, &
    wR, wI
  integer :: n, lastY, lastZ
      
  call MPI_Init_thread(MPI_THREAD_FUNNELED,provided,ierr)
  
  call GETARG(1,inputfileLS)
  call GETARG(2,inputfileSS)
  call GETARG(3,inputfileGM)
  call GetArguments(nthreads)
  
  call largeScales%init(inputfileLS, .true.) 
  call smallScales%init(inputfileSS, .false.)

  call enrich%init(smallScales,largeScales,inputfileGM)

  enrich%uhatR = real(nrank,rkind)
  enrich%uhatI = 0.d0
  enrich%vhatR = 0.d0
  enrich%vhatI = 0.d0
  enrich%whatR = 0.d0
  enrich%whatI = 0.d0

  enrich%kx = 0.d0
  enrich%ky = 0.d0
  enrich%kz = 0.d0

  if (allocated(enrich%haloBuffY)) deallocate(enrich%haloBuffY)
  if (allocated(enrich%haloBuffZ)) deallocate(enrich%haloBuffZ)
  call enrich%sendRecvHaloModes(enrich%modeData, 'y', enrich%haloBuffY, lastY)
  call enrich%sendRecvHaloModes(enrich%modeData, 'z', enrich%haloBuffZ, &
    lastZ, enrich%haloBuffY)

  call assert(size(enrich%haloBuffY,1) == lastY,&
    'size(enrich%haloBuffY,1) == lastY')
  call assert(size(enrich%haloBuffZ,1) == lastZ,&
    'size(enrich%haloBuffZ,1) == lastZ')

  allocate(x(size(enrich%haloBuffY,1) + size(enrich%haloBuffZ,1)))
  allocate(y(size(enrich%haloBuffY,1) + size(enrich%haloBuffZ,1)))
  allocate(z(size(enrich%haloBuffY,1) + size(enrich%haloBuffZ,1)))

  allocate(kx(size(enrich%haloBuffY,1) + size(enrich%haloBuffZ,1)))
  allocate(ky(size(enrich%haloBuffY,1) + size(enrich%haloBuffZ,1)))
  allocate(kz(size(enrich%haloBuffY,1) + size(enrich%haloBuffZ,1)))

  allocate(uR(size(enrich%haloBuffY,1) + size(enrich%haloBuffZ,1)))
  allocate(uI(size(enrich%haloBuffY,1) + size(enrich%haloBuffZ,1)))
  allocate(vR(size(enrich%haloBuffY,1) + size(enrich%haloBuffZ,1)))
  allocate(vI(size(enrich%haloBuffY,1) + size(enrich%haloBuffZ,1)))
  allocate(wR(size(enrich%haloBuffY,1) + size(enrich%haloBuffZ,1)))
  allocate(wI(size(enrich%haloBuffY,1) + size(enrich%haloBuffZ,1)))

  x = [enrich%haloBuffY(:,1), enrich%haloBuffZ(:,1)]
  y = [enrich%haloBuffY(:,2), enrich%haloBuffZ(:,2)]
  z = [enrich%haloBuffY(:,3), enrich%haloBuffZ(:,3)]

  kx = [enrich%haloBuffY(:,4), enrich%haloBuffZ(:,4)]
  ky = [enrich%haloBuffY(:,5), enrich%haloBuffZ(:,5)]
  kz = [enrich%haloBuffY(:,6), enrich%haloBuffZ(:,6)]
  
  uR = [enrich%haloBuffY(:,7), enrich%haloBuffZ(:,7)]
  uI = [enrich%haloBuffY(:,8), enrich%haloBuffZ(:,8)]
  vR = [enrich%haloBuffY(:,9), enrich%haloBuffZ(:,9)]
  vI = [enrich%haloBuffY(:,10), enrich%haloBuffZ(:,10)]
  wR = [enrich%haloBuffY(:,11), enrich%haloBuffZ(:,11)]
  wI = [enrich%haloBuffY(:,12), enrich%haloBuffZ(:,12)]

  call enrich%dumpData(x,y,z,kx,ky,kz,uR,uI,vR,vI,wR,wI)

  do n = 1,nproc
    if (nrank == n-1) then
      print*, "rank", nrank
      print*, "(ymin,zmin):", &
        enrich%QHgrid%yE(1), enrich%QHgrid%zE(1) 
      print*, "(ymax,zmax):", &
        enrich%QHgrid%yE(enrich%QHgrid%gpC%xsz(2)+1), &
        enrich%QHgrid%zE(enrich%QHgrid%gpC%xsz(3)+1) 
      print*, "----------------------------------"
      print*, " "
    end if
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
  end do

  deallocate(x,y,z,kx,ky,kz,uR,uI,vR,vI,wR,wI)
  call enrich%destroy()
  call largeScales%destroy()
  call smallScales%destroy()

  call MPI_Finalize(ierr)
end program 
