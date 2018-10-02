! Template for PadeOps

#include "pblwt_files/initialize.F90"       
#include "pblwt_files/temporalHook.F90"  

program pblwt
    use mpi
    use kind_parameters,  only: clen
    use IncompressibleGrid, only: igrid
    !use IncompressibleGridWallM, only: igridWallM
    use temporalhook, only: doTemporalStuff
    use timer, only: tic, toc
    use exits, only: message

    implicit none

    !type(igridWallM), allocatable, target :: igp
    type(igrid), allocatable, target :: igp
    character(len=clen) :: inputfile
    integer :: ierr

    call MPI_Init(ierr)               !<-- Begin MPI

    call GETARG(1,inputfile)          !<-- Get the location of the input file

    allocate(igp)                     !<-- Initialize hit_grid with defaults

    call igp%init(inputfile)          !<-- Properly initialize the hit_grid solver (see hit_grid.F90)
  
    call igp%start_io(.true.)                !<-- Start I/O by creating a header file (see io.F90)

    call igp%printDivergence()
  
    call tic() 
    do while ((igp%tsim < igp%tstop) .and. (igp%step < igp%nsteps))
       
       call igp%timeAdvance()     !<-- Time stepping scheme + Pressure Proj. (see igridWallM.F90 or igrid.F90)
       call doTemporalStuff(igp)     !<-- Go to the temporal hook (see temporalHook.F90)
       
    end do 
 
    call igp%finalize_io()                  !<-- Close the header file (wrap up i/o)

    call igp%destroy_grid()                !<-- Destroy the IGRID derived type 
   

    deallocate(igp)                   !<-- Deallocate all the memory associated with scalar defaults
    
    call MPI_Finalize(ierr)           !<-- Terminate MPI 

end program
