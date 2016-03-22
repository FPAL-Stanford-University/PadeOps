#include "PlasticRelaxation_files/hooks.F90"

program PlasticRelaxation

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

    ! Time advance
    call sgp%simulate()
        
    ! Destroy everythin before ending
    call sgp%destroy()

    ! End the run
    call MPI_Finalize(ierr)

end program
