! Template for PadeOps
! Grid - igrid
! Problem - HIT

#include "hit_files/meshgen.F90"       
#include "hit_files/io.F90"            
#include "hit_files/initfields.F90"    
#include "hit_files/temporalHook.F90"  

program hitcd
    use mpi
    use kind_parameters,  only: rkind,clen,stdout,stderr
    use IncompressibleGrid, only: igrid
    use hitCD_IO, only: start_io, finalize_io
    use constants, only: half 
    use temporalhook, only: doTemporalStuff
    use timer, only: tic, toc 
    implicit none

    type(igrid), allocatable, target :: igp
    character(len=clen) :: inputfile
    integer :: ierr

    call MPI_Init(ierr)               !<-- Begin MPI

    call GETARG(1,inputfile)          !<-- Get the location of the input file

    allocate(igp)                     !<-- Initialize igrid with defaults

    call igp%init(inputfile)          !<-- Properly initialize the igrid solver (see igrid.F90)
    
    call start_io(igp)                !<-- Start I/O by creating a header file (see io.F90)


    call igp%printDivergence()
    do while (igp%tsim < igp%tstop) 
        
        call igp%AdamsBashforth()     !<-- Time stepping scheme + Pressure Proj. (see igrid.F90)
        
        call doTemporalStuff(igp)     !<-- Go to the temporal hook (see temporalHook.F90)
        
        igp%tsim = igp%tsim + igp%dt
        igp%step = igp%step + 1
    end do 
    
    call finalize_io                  !<-- Close the header file (wrap up i/o)

    call igp%destroy()                !<-- Destroy the IGRID derived type 
    
    deallocate(igp)                   !<-- Deallocate all the memory associated with scalar defaults
    
    call MPI_Finalize(ierr)           !<-- Terminate MPI 

end program
