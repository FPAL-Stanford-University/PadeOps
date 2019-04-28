! Template for PadeOps

#include "EkmanLayerDNS_files/io.F90"       
#include "EkmanLayerDNS_files/initialize.F90"       
#include "EkmanLayerDNS_files/temporalHook.F90"  

program EkmanLayerDNS
    use mpi
    use kind_parameters,  only: clen
    use IncompressibleGrid, only: igrid
    use temporalhook, only: doTemporalStuff, initialize_processing, finalize_processing
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

    call tic()

    call initialize_processing(igp%nz, inputfile)
    do while (igp%tsim < igp%tstop) 
       
       call igp%timeAdvance()     !<-- Time stepping scheme + Pressure Proj. (see igrid.F90)
       call doTemporalStuff(igp)     !<-- Go to the temporal hook (see temporalHook.F90)
    end do 

    call finalize_processing()
    call igp%finalize_io()                  !<-- Close the header file (wrap up i/o)

    call igp%destroy()                !<-- Destroy the IGRID derived type 
   

    deallocate(igp)                   !<-- Deallocate all the memory associated with scalar defaults
    
    call MPI_Finalize(ierr)           !<-- Terminate MPI 

end program
