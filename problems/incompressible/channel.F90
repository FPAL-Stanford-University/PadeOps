! Template for PadeOps
! Grid - channel_grid
! Problem - HIT

#include "channel_files/initialize.F90"       
#include "channel_files/io.F90"            
#include "channel_files/temporalHook.F90"  

program EkmanBL
    use mpi
    use kind_parameters,  only: rkind,clen,stdout,stderr
    use IncompressibleGridNP, only: igrid
    use EkmanBL_IO, only: start_io, finalize_io
    use constants, only: half 
    use temporalhook, only: doTemporalStuff
    use timer, only: tic, toc
    use allStatistics, only: init_stats, finalize_stats

    implicit none

    type(igrid), allocatable, target :: igp
    character(len=clen) :: inputfile
    integer :: ierr

    call MPI_Init(ierr)               !<-- Begin MPI

    call GETARG(1,inputfile)          !<-- Get the location of the input file

    allocate(igp)                     !<-- Initialize hit_grid with defaults

    call igp%init(inputfile)          !<-- Properly initialize the hit_grid solver (see hit_grid.F90)
  
    call start_io(igp)                !<-- Start I/O by creating a header file (see io.F90)

    call igp%printDivergence()
  
    call init_stats(igp)  

    call tic() 
    do while (igp%tsim < igp%tstop) 
       
       call igp%AdamsBashforth()     !<-- Time stepping scheme + Pressure Proj. (see hit_grid.F90)
      
       igp%step = igp%step + 1 
       igp%tsim = igp%tsim + igp%dt
       call doTemporalStuff(igp)     !<-- Go to the temporal hook (see temporalHook.F90)
       
    end do 
    
    call finalize_io                  !<-- Close the header file (wrap up i/o)

    call igp%destroy()                !<-- Destroy the IGRID derived type 
   
    call finalize_stats()

    deallocate(igp)                   !<-- Deallocate all the memory associated with scalar defaults
    
    call MPI_Finalize(ierr)           !<-- Terminate MPI 

end program
