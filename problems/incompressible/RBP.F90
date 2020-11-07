! Template for PadeOps

#include "RBP_files/initialize.F90"       
#include "RBP_files/temporalHook.F90"  

program RBPGrowth
    use mpi
    use kind_parameters,  only: clen, rkind
    use IncompressibleGrid, only: igrid
    use temporalhook, only: doTemporalStuff
    use timer, only: tic, toc
    use exits, only: message
    use RBP_parameters
    implicit none

    type(igrid), allocatable, target :: igp
    character(len=clen) :: inputfile
    integer :: ierr

    call MPI_Init(ierr)               !<-- Begin MPI

    call GETARG(1,inputfile)          !<-- Get the location of the input file

    allocate(igp)                     !<-- Initialize hit_grid with defaults

    call igp%init(inputfile)          !<-- Properly initialize the hit_grid solver (see hit_grid.F90)
    call SetTemperatureBC_RBP(igp,inputfile)

    call igp%start_io(.true.)         !<-- Start I/O by creating a header file (see io.F90)

    call igp%printDivergence()
   
    ! Allocate memory for fringe 
    call igp%fringe_x%allocateTargetArray_Cells(utarget0)                
    call igp%fringe_x%allocateTargetArray_Cells(vtarget0)                
    call igp%fringe_x%allocateTargetArray_Edges(wtarget0)                
    call igp%fringe_x%allocateTargetArray_Cells(Ttarget0)                
 
    ! Base state at fringe
    call setBaseState(igp%mesh(:,:,:,3),utarget0,vtarget0,wtarget0,Ttarget0,inputfile)

    ! Link Fringe Targets
    call igp%fringe_x%associateFringeTargets(utarget0, vtarget0, wtarget0) 
    call igp%fringe_x%associateFringeTarget_scalar(Ttarget0)

    call tic() 
    do while (igp%tsim < igp%tstop) 
       
       call igp%timeAdvance()        !<-- Time stepping scheme + Pressure Proj. (see igrid.F90)
       call doTemporalStuff(igp)     !<-- Go to the temporal hook (see temporalHook.F90)
       
    end do 
 
    call igp%finalize_io()            !<-- Close the header file (wrap up i/o)

    call igp%destroy()                !<-- Destroy the IGRID derived type 

    deallocate(igp)                   !<-- Deallocate all the memory associated with scalar defaults
    
    call MPI_Finalize(ierr)           !<-- Terminate MPI 

end program
