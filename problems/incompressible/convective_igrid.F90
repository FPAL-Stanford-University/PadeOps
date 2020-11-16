! Template for PadeOps

#include "convective_igrid_files/initialize.F90"       
#include "convective_igrid_files/temporalHook.F90"  

program convective_igrid
    use mpi
    use kind_parameters,  only: clen
    use IncompressibleGrid, only: igrid
    use temporalhook, only: doTemporalStuff
    use timer, only: tic, toc
    use exits, only: message
    use budgets_time_avg_mod, only: budgets_time_avg  
    use budgets_xy_avg_mod,   only: budgets_xy_avg  

    implicit none

    type(budgets_time_avg) :: budg_tavg
    type(budgets_xy_avg)   :: budg_xyavg
    type(igrid), allocatable, target :: igp
    character(len=clen) :: inputfile
    integer :: ierr

    call MPI_Init(ierr)               !<-- Begin MPI

    call GETARG(1,inputfile)          !<-- Get the location of the input file

    allocate(igp)                     !<-- Initialize hit_grid with defaults

    call igp%init(inputfile)          !<-- Properly initialize the hit_grid solver (see hit_grid.F90)
  
    call igp%start_io(.false.)                !<-- Start I/O by creating a header file (see io.F90)
    
    call igp%printDivergence()
  
    call budg_xyavg%init(inputfile, igp)   !<-- Budget class initialization 
    call budg_tavg%init( inputfile, igp)   !<-- Budget class initialization 

    call tic() 
    do while (igp%tsim < igp%tstop) 
       
       call igp%timeAdvance()     !<-- Time stepping scheme + Pressure Proj. (see igridWallM.F90)
       call budg_xyavg%doBudgets()       
       call budg_tavg%doBudgets()       
       call doTemporalStuff(igp)     !<-- Go to the temporal hook (see temporalHook.F90)
       
    end do 
 
    call budg_xyavg%doBudgets(.true.)      !<-- Force dump budget information if active
    call budg_tavg%doBudgets(.true.)       !<-- Force dump budget information if active

    call igp%finalize_io()                  !<-- Close the header file (wrap up i/o)

    call budg_xyavg%destroy()          !<-- release memory taken by the budget class 
    call budg_tavg%destroy()           !<-- release memory taken by the budget class 

    call igp%destroy()                !<-- Destroy the IGRID derived type 
   
    deallocate(igp)                   !<-- Deallocate all the memory associated with scalar defaults
    
    call MPI_Finalize(ierr)           !<-- Terminate MPI 

end program
