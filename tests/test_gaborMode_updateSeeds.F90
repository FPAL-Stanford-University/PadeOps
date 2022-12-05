#include "../problems/gabor/hitEnrich_files/initialize.F90"
program test_gaborMode_renderVelMPI
  use kind_parameters,       only: clen, rkind
  use enrichmentMod,         only: enrichmentOperator, nthreads, xDom, yDom, zDom
  use IncompressibleGrid,    only: igrid
  use auxiliary_openmp_subs, only: GetArguments
  use mpi
  use exits,                 only: message
  use decomp_2d,             only: nrank
  use GaborModeRoutines,     only: initializeSeeds, getIdx, updateSeeds
  implicit none

  character(len=clen) :: inputfileLS, inputfileSS, inputfileGM 
  character(len=clen) :: datadir, fname, outputdir
  integer :: ioUnit, ierr, provided
  type(igrid) :: largeScales, smallScales
  type(enrichmentOperator) :: enrich
  real(rkind) :: Lx, Ly, Lz, xWindow
  integer :: idxOld, idxCurrent
  real(rkind), dimension(7) :: seed
  integer :: nx, ny
  integer :: ist, jst, kst, i, j, k
      
  call MPI_Init_thread(MPI_THREAD_FUNNELED,provided,ierr)
  
  call GETARG(1,inputfileLS)
  call GETARG(2,inputfileSS)
  call GETARG(3,inputfileGM)
  call GetArguments(nthreads,4)
  
  call largeScales%init(inputfileLS, .true.) 
  call smallScales%init(inputfileSS, .false.)

  call enrich%init(smallScales,largeScales,inputfileGM)

  ist = enrich%QHgrid%gpC%xst(1)
  jst = enrich%QHgrid%gpC%xst(2)
  kst = enrich%QHgrid%gpC%xst(3) 

  nx = enrich%QHgrid%gpC%xsz(1)
  ny = enrich%QHgrid%gpC%ysz(2)

  call initializeSeeds(seed,ist,jst,kst,nx,ny)

  idxOld = getIdx(ist,jst,kst,nx,ny) - 1
  do k = 1,enrich%QHgrid%gpC%xsz(3)
    do j = 1,enrich%QHgrid%gpC%xsz(2)
      do i = 1,enrich%QHgrid%gpC%xsz(1)
        idxCurrent = getIdx(ist + i - 1, jst + j - 1, kst + k - 1, nx, ny)
        call updateSeeds(seed,idxOld,idxCurrent)
        idxOld = idxCurrent
        print*, ist + i - 1, jst + j - 1, kst + k - 1, nint(seed(1))
      end do
    end do
  end do

  call enrich%destroy()
  call largeScales%destroy()
  call smallScales%destroy()

  call MPI_Finalize(ierr)
end program 
