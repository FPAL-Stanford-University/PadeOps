#include "test_pTequilibrium_files/hooks.F90"

program test_ptequilibrium

    use kind_parameters,  only: clen
    use SolidGrid,        only: sgrid
    implicit none

    type(sgrid) :: sgp
    character(len=clen) :: inputfile
    integer :: ierr 

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

    write(*,*) 'Mat 1 P, T before: ', sgp%mix%material(1)%p(1,1,1), sgp%mix%material(1)%T(1,1,1) 
    write(*,*) 'Mat 2 P, T before: ', sgp%mix%material(2)%p(1,1,1), sgp%mix%material(2)%T(1,1,1) 

    
    sgp%mix%material(1)%p(1,1,1) = 1.0e-1+sgp%mix%material(1)%p(1,1,1) 
    sgp%mix%material(1)%T(1,1,1) = 2.0*sgp%mix%material(1)%T(1,1,1) 

    call sgp%mix%equilibratePressureTemperature(sgp%rho, sgp%e, sgp%p, sgp%T)

    write(*,*) 'Mat 1 P, T after: ', sgp%mix%material(1)%p(1,1,1), sgp%mix%material(1)%T(1,1,1) 
    write(*,*) 'Mat 2 P, T after: ', sgp%mix%material(2)%p(1,1,1), sgp%mix%material(2)%T(1,1,1) 
    
    ! Destroy everythin before ending
    call sgp%destroy()
    write(*,*) 'Done destroy'

    ! End the run
    call MPI_Finalize(ierr)
    write(*,*) 'Done finalize'

end program
