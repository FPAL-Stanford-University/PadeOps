#include "ManufacturedSolution_files/hooks.F90"

program ManufacturedSolution

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
    print *, 'A0'   
    
    ! Initialize the grid object
    call sgp%init(inputfile)
    print *, 'A'   

    ! Time advance
    call sgp%simulate()
    print *, 'B'   
     
    ! Destroy everythin before ending
    call sgp%destroy()
    print *, 'C'   

    ! End the run
    call MPI_Finalize(ierr)

end program
