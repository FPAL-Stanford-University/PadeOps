#include "CavityCollapse1D_files/hooks.F90"

program CavityCollapse1D

    use kind_parameters,  only: clen
    use SolidGrid,        only: sgrid
    use decomp_2d,        only: nrank
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
    if(nrank==0) write (*,*) 'Done init'

    ! Time advance
    call sgp%simulate()
    if(nrank==0) write (*,*) 'Done simulate'
        
    ! Destroy everythin before ending
    call sgp%destroy()

    ! End the run
    call MPI_Finalize(ierr)

end program
