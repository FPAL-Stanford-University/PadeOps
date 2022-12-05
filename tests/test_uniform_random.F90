module update_seed_mod
  implicit none
  contains
    subroutine update_seed(seed,iter)
      integer, intent(inout) :: seed
      integer, intent(in) :: iter
      seed = 1294050*iter
    end subroutine
end module

program test_uniform_random
  use mpi
  use kind_parameters, only: rkind
  use random, only: uniform_random
  use update_seed_mod, only: update_seed
  implicit none

  real(rkind), dimension(4) :: rand1
  integer :: seed1, seed2, n, niter
  integer :: ierr
  character(len=3) :: niterstr

  call MPI_Init(ierr)

  call GETARG(1,niterstr)
  
  seed1 = 97059
  seed2 = 48538
  read (niterstr,'(I3)') niter
  do n = 1,niter
    call uniform_random(rand1,0.d0,1.d0,seed1)
    call uniform_random(rand2,0.d0,1.d0,seed2)
    print*, "rand1:", rand1
    print*, "rand2:", rand2
    call update_seed(seed1,n)
  end do
  call MPI_Finalize(ierr)
end program test_uniform_random
