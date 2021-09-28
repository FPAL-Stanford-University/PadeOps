! Template for PadeOps

#include "AD_Coriolis_files/initialize.F90"       
#include "AD_Coriolis_files/temporalHook.F90"  

program AD_Coriolis
    use mpi
    use kind_parameters,  only: clen, rkind
    use IncompressibleGrid, only: igrid
    use temporalhook, only: doTemporalStuff
    use timer, only: tic, toc
    use exits, only: message
    use AD_Coriolis_parameters 
    use budgets_time_avg_mod, only: budgets_time_avg  

    implicit none

    type(igrid), allocatable, target :: igp
    character(len=clen) :: inputfile
    integer :: ierr
    type(budgets_time_avg) :: budg_tavg

    call MPI_Init(ierr)                                                 !<-- Begin MPI

    call GETARG(1,inputfile)                                            !<-- Get the location of the input file

    allocate(igp)                                                       !<-- Initialize hit_grid with defaults

    call igp%init(inputfile)                                            !<-- Properly initialize the hit_grid solver (see hit_grid.F90)
    call igp%start_io(.true.)                                           !<-- Start I/O by creating a header file (see io.F90)
    call igp%printDivergence()
 
    ! Budgets
    call budg_tavg%init(inputfile, igp)   !<-- Budget class initialization 
    
    ! Fringe associations for non-periodic BCs in x
    call igp%fringe_x%allocateTargetArray_Cells(utarget)                !<-- Allocate target array of appropriate size
    call igp%fringe_x%allocateTargetArray_Cells(vtarget)                !<-- Allocate target array of appropriate size
    call igp%fringe_x%allocateTargetArray_Edges(wtarget)                !<-- Allocate target array of appropriate size
    call init_fringe_targets(inputfile, igp%mesh)  
    call igp%fringe_x%associateFringeTargets(utarget, vtarget, wtarget) !<-- Link the target velocity array to igp 

    ! NOTE: Beyond this point, the target arrays can change in time, In case of
    ! concurrent simulations where inflow is turbulent, the utarget, vtarget and
    ! wtarget are simply set equal to the u, v, and w arrays of the concurrent
    ! simulation (at the end of each time step).  

    call tic() 
    do while (igp%tsim < igp%tstop) 
       
       call igp%timeAdvance()                                           !<-- Time stepping scheme + Pressure Proj. (see igrid.F90)
       call budg_tavg%doBudgets()       
       call doTemporalStuff(igp)                                        !<-- Go to the temporal hook (see temporalHook.F90)
       
    end do 
    
    call budg_tavg%destroy()           !<-- release memory taken by the budget class 
 
    call igp%finalize_io()                                              !<-- Close the header file (wrap up i/o)

    call igp%destroy()                                                  !<-- Destroy the IGRID derived type 
   

    deallocate(igp)                                                     !<-- Deallocate all the memory associated with scalar defaults
    
    call MPI_Finalize(ierr)                                             !<-- Terminate MPI 

end program
