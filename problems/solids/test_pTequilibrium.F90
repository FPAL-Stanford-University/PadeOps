#include "test_pTequilibrium_files/hooks.F90"

program test_ptequilibrium

    use kind_parameters,  only: clen
    use SolidGrid,        only: sgrid
    implicit none

    type(sgrid) :: sgp
    character(len=clen) :: inputfile
    integer :: ierr, Nx 

    ! Start MPI
    call MPI_Init(ierr)

    ! Get file location 
    call GETARG(1,inputfile)
    
    ! Initialize the grid object
    call sgp%init(inputfile)
    write(*,*) 'Done init'

    !! Time advance
    !call sgp%simulate()
    !write(*,*) 'Done simulate'

    Nx = 65*2
    write(*,*) 'Vol fractions', sgp%mix%material(1)%VF(63:67,1,1)
    write(*,*) 'VF, Ys at location: ', sgp%mix%material(1)%VF(Nx/2,1,1), sgp%mix%material(1)%Ys(Nx/2,1,1) 

    write(*,*) 'Mat 1 P, T   before: ', sgp%mix%material(1)%p(Nx/2,1,1), sgp%mix%material(1)%T(Nx/2,1,1) 
    write(*,*) 'Mat 2 P, T   before: ', sgp%mix%material(2)%p(Nx/2,1,1), sgp%mix%material(2)%T(Nx/2,1,1) 
    write(*,*) 'Mat 1 eh vf  before: ', sgp%mix%material(1)%energy(Nx/2,1,1), sgp%mix%material(1)%VF(Nx/2,1,1) 
    write(*,*) 'Mat 2 eh vf  before: ', sgp%mix%material(2)%energy(Nx/2,1,1), sgp%mix%material(2)%VF(Nx/2,1,1) 

    write(*,*) 'Mix e, rho before: ', sgp%e(Nx/2,1,1), sgp%rho(Nx/2,1,1) 

    !! perturbation 1
    !sgp%mix%material(1)%p(Nx/2,1,1) = 0.0e-1+sgp%mix%material(1)%p(Nx/2,1,1) 
    !sgp%mix%material(1)%T(Nx/2,1,1) = 2.0*sgp%mix%material(1)%T(Nx/2,1,1) 
    !call sgp%post_bc()
    !write(*,*) 'Perturbed Mat 1 p: ', sgp%mix%material(1)%p(Nx/2,1,1)
    !write(*,*) 'Perturbed Mat 1 T: ', sgp%mix%material(1)%T(Nx/2,1,1)

    ! perturbation 2
    sgp%mix%material(1)%energy(:,1,1) = 2.0*sgp%mix%material(1)%energy(:,1,1) 
    sgp%mix%material(2)%energy(:,1,1) = 2.0*sgp%mix%material(2)%energy(:,1,1) 
    write(*,*) 'Perturbed Mat 1 eh: ', sgp%mix%material(1)%energy(Nx/2,1,1)
    write(*,*) 'Perturbed Mat 2 eh: ', sgp%mix%material(2)%energy(Nx/2,1,1)

    call sgp%mix%equilibratePressureTemperature(sgp%rho, sgp%e)

    write(*,*) 'Mat 1 P, T   after: ', sgp%mix%material(1)%p(Nx/2,1,1), sgp%mix%material(1)%T(Nx/2,1,1) 
    write(*,*) 'Mat 2 P, T   after: ', sgp%mix%material(2)%p(Nx/2,1,1), sgp%mix%material(2)%T(Nx/2,1,1) 
    write(*,*) 'Mat 1 eh vf  after: ', sgp%mix%material(1)%energy(Nx/2,1,1), sgp%mix%material(1)%VF(Nx/2,1,1) 
    write(*,*) 'Mat 2 eh vf  after: ', sgp%mix%material(2)%energy(Nx/2,1,1), sgp%mix%material(2)%VF(Nx/2,1,1) 
    
    ! Destroy everythin before ending
    call sgp%destroy()
    write(*,*) 'Done destroy'

    ! End the run
    call MPI_Finalize(ierr)
    write(*,*) 'Done finalize'

end program
