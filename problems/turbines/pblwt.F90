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
    use budgets_time_avg_mod, only: budgets_time_avg  
    use budgets_xy_avg_mod,   only: budgets_xy_avg  

    implicit none

    !type(igridWallM), allocatable, target :: igp
    type(igrid), allocatable, target :: igp
    character(len=clen) :: inputfile
    integer :: ierr
    type(budgets_time_avg) :: budg_tavg
    type(budgets_xy_avg)   :: budg_xyavg1, budg_xyavg2

    call MPI_Init(ierr)               !<-- Begin MPI

    call GETARG(1,inputfile)          !<-- Get the location of the input file

    allocate(igp)                     !<-- Initialize hit_grid with defaults

    call igp%init(inputfile)          !<-- Properly initialize the hit_grid solver (see hit_grid.F90)
  
    call igp%start_io(.true.)                !<-- Start I/O by creating a header file (see io.F90)

    call igp%printDivergence()
  
    call budg_xyavg1%init('input_budgxy1.dat', igp)   !<-- Budget class initialization 
    call budg_xyavg2%init('input_budgxy2.dat', igp)   !<-- Budget class initialization 
  
    call budg_tavg%init(inputfile, igp)   !<-- Budget class initialization 
    
    call tic() 
    do while ((igp%tsim < igp%tstop) .and. (igp%step < igp%nsteps))
       
       call igp%timeAdvance()     !<-- Time stepping scheme + Pressure Proj. (see igridWallM.F90 or igrid.F90)
       call budg_xyavg1%doBudgets()       
       call budg_xyavg2%doBudgets()       
       call budg_tavg%doBudgets()       
       call doTemporalStuff(igp)     !<-- Go to the temporal hook (see temporalHook.F90)
       
    end do 
 
    call budg_xyavg1%doBudgets(.true.)      !<-- Force dump budget information if active
    call budg_xyavg2%doBudgets(.true.)      !<-- Force dump budget information if active

    call budg_tavg%doBudgets(.true.)       !<-- Force dump budget information if active

    call igp%finalize_io()                  !<-- Close the header file (wrap up i/o)
    
    call budg_xyavg1%destroy()          !<-- release memory taken by the budget class 
    call budg_xyavg2%destroy()          !<-- release memory taken by the budget class 

    call budg_tavg%destroy()           !<-- release memory taken by the budget class 

    call igp%destroy()                !<-- Destroy the IGRID derived type 
   

    deallocate(igp)                   !<-- Deallocate all the memory associated with scalar defaults
    
    call MPI_Finalize(ierr)           !<-- Terminate MPI 

end program
