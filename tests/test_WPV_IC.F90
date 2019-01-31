! Template for PadeOps

#include "test_WPV_IC_files/initialize.F90"       
#include "test_WPV_IC_files/io.F90"

program test_WPV_IC_igrid
    use mpi
    use kind_parameters,  only: clen
    use IncompressibleGrid, only: igrid
    use timer, only: tic, toc
    use exits, only: message

    implicit none

    type(igrid), allocatable, target :: igp
    character(len=clen) :: inputfile
    integer :: ierr

    call MPI_Init(ierr)               !<-- Begin MPI

    call GETARG(1,inputfile)          !<-- Get the location of the input file

    allocate(igp)                     !<-- Initialize hit_grid with defaults

    call igp%init(inputfile)          !<-- Properly initialize the hit_grid solver (see hit_grid.F90)
  
    call igp%start_io(.true.)                !<-- Start I/O by creating a header file (see io.F90)
    
    call igp%printDivergence()
  
    call igp%finalize_io()                  !<-- Close the header file (wrap up i/o)

    call igp%destroy()                !<-- Destroy the IGRID derived type 
   

    deallocate(igp)                   !<-- Deallocate all the memory associated with scalar defaults
    
    call MPI_Finalize(ierr)           !<-- Terminate MPI 

end program
