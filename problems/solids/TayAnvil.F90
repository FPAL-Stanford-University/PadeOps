#include "TayAnvil_files/hooks.F90"

program TayAnvil

    use kind_parameters,  only: clen
    use SolidGrid,        only: sgrid
    use decomp_2d,        only: nrank
    use mpi
    implicit none

    type(sgrid) :: sgp
    character(len=clen) :: inputfile
    integer :: ierr 

    ! Start MPI
    call MPI_Init(ierr)
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    if(nrank==0) write (*,*) 'Done mpi init'
    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    ! Get file location 
    call GETARG(1,inputfile)
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    if(nrank==0) write (*,*) 'Done getarg'
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    
    ! Initialize the grid object
    call sgp%init(inputfile)
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    if(nrank==0) write (*,*) 'Done init'
    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    ! Time advance
    call sgp%simulate()
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    if(nrank==0) write (*,*) 'Done simulate'
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
        
    ! Destroy everythin before ending
    call sgp%destroy()

    ! End the run
    call MPI_Finalize(ierr)

end program
