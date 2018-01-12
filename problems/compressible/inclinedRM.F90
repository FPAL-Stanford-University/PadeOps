#include "inclinedRM_files/hooks.F90"

program inclinedRM

    use kind_parameters,  only: rkind,clen
    use CompressibleGrid, only: cgrid
    implicit none

    type(cgrid) :: cgp
    character(len=clen) :: inputfile
    integer :: ierr 

    ! Start MPI
    call MPI_Init(ierr)

    ! Get file location 
    call GETARG(1,inputfile)
    
    ! Initialize the grid object
    call cgp%init(inputfile)

    ! Time advance
    call cgp%simulate()
        
    ! Destroy everythin before ending
    call cgp%destroy()

    ! End the run
    call MPI_Finalize(ierr)

end program
