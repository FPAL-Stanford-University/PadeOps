#include "../problems/gabor/hitEnrich_files/initialize.F90"
module testMod
  use mpi
  use kind_parameters, only: rkind
  use decomp_2d, only: nrank, nproc
  use enrichmentMod, only: enrichmentOperator
  implicit none

  contains
    subroutine printNonZeroModeInfo(eop)
      type(enrichmentOperator), intent(in) :: eop
      integer :: m, n, ierr, id
      integer, dimension(:), allocatable :: idx
      real(rkind) :: tol = 1.d-13

      id = 0
      do m = 1,nproc
          if (m-1 == nrank) then
              print*, "PE", nrank, "B"
              do n = 1,size(eop%uhatR)
                  if (abs(eop%uhatR(n)) > tol) id = id + 1
              end do
              allocate(idx(id))
              id = 0
              do n = 1,size(eop%uhatR)
                  if (abs(eop%uhatR(n)) > tol) then
                    id = id + 1
                    idx(id) = n
                  end if
              end do
              do id = 1,size(idx)
                  do n = 1,12
                      print*, eop%modeData(idx(id),n)
                  end do
                  print*, "-------------------------------------"
              end do
              deallocate(idx)
          end if
          call MPI_Barrier(MPI_COMM_WORLD,ierr)
      end do
    end subroutine
    subroutine zeroOnly(eop,idx)
      type(enrichmentOperator), intent(inout) :: eop
      integer, intent(in) :: idx
      integer :: n
      
      do n = 1,size(eop%modeData,1)
        if (n == idx) then
          eop%modeData(n,7:12) =  0.d0
        else
          continue
        end if
      end do 
    end subroutine
    subroutine zeroAllBut(eop,idx1,idx2)
      type(enrichmentOperator), intent(inout) :: eop
      integer, intent(in) :: idx1
      integer, intent(in), optional :: idx2
      integer :: n
      if (present(idx2)) then
        do n = 1,size(eop%modeData,1)
          if (n >= idx1 .and. n <= idx2) then 
            continue
          else
            eop%modeData(n,7:12) = 0.d0
          end if
        end do
      else
        do n = 1,size(eop%modeData,1)
          if (n == idx1) then 
            continue
          else
            eop%modeData(n,7:12) = 0.d0
          end if
        end do
      end if
    end subroutine
    subroutine printModeInfo(eop,idx)
      type(enrichmentOperator), intent(in) :: eop
      integer, intent(in) :: idx
      integer :: m, n, ierr
      do m = 1,nproc
        if (m-1 == nrank) then
          print*, "PE", nrank, "A"
          do n = 1,12
            print*, eop%modeData(idx,n)
          end do
          print*, "-------------------------------------"
        else 
          continue
        end if
        call MPI_Barrier(MPI_COMM_WORLD,ierr)
      end do
    end subroutine
end module

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
  use testMod
  implicit none

  character(len=clen) :: inputfileLS, inputfileSS, inputfileGM 
  integer :: ierr, provided
  type(igrid) :: largeScales, smallScales
  type(enrichmentOperator) :: enrich
  real(rkind), dimension(:,:,:), allocatable :: u1, v1, w1
  real(rkind) :: tol = 1.d-13
  real(rkind) :: Lz
  integer :: n, m, id, idx, idx1, idx2
  character(len=2) :: idx1str, idx2str
      
  call MPI_Init_thread(MPI_THREAD_FUNNELED,provided,ierr)
  
  call GETARG(1,inputfileLS)
  call GETARG(2,inputfileSS)
  call GETARG(3,inputfileGM)
  call GETARG(4,idx1str)
  call GETARG(5,idx2str)
  call GetArguments(nthreads,6)

  read(idx1str,'(I2)') idx1
  read(idx2str,'(I2)') idx2
  
  call largeScales%init(inputfileLS, .true.) 
  call smallScales%init(inputfileSS, .false.)
  call enrich%init(smallScales,largeScales,inputfileGM)
  Lz = zDom(2) - zDom(1)

  ! Set uniform flow in one coordinate direction
  largeScales%v       = two*pi !enrich%smallScales%dy/enrich%dt 
  largeScales%u       = 0.d0
  largeScales%wC      = 0.d0
  largeScales%duidxjC = 0.d0
  
  call zeroAllBut(enrich,idx1,idx2)
  do n = idx1,idx2
    call printModeInfo(enrich,n)
  end do
  
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

  call printNonZeroModeInfo(enrich)
  
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
