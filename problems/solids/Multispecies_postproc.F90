#include "Multispecies_postproc_files/hooks.F90"

program Multispecies_postproc

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

    ! postprocess
    call sgp%postproc()
    if(nrank==0) write (*,*) 'Done postproc'
        
    ! Destroy everythin before ending
    call sgp%destroy()

    ! End the run
    call MPI_Finalize(ierr)

end program
