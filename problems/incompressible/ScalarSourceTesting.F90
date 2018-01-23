! Template for PadeOps

#include "ScalarSourceTesting_files/initialize.F90"       
#include "ScalarSourceTesting_files/temporalHook.F90"  

program ScalarSourceTesting
    use mpi
    use kind_parameters,  only: rkind, clen
    use IncompressibleGrid, only: igrid
    use temporalhook, only: doTemporalStuff
    use timer, only: tic, toc
    use exits, only: message

    implicit none

    type(igrid), allocatable, target :: igp
    character(len=clen) :: inputfile
    integer :: ierr
    real(rkind), dimension(:,:,:), allocatable :: utarget, vtarget, wtarget

    call MPI_Init(ierr)               !<-- Begin MPI

    call GETARG(1,inputfile)          !<-- Get the location of the input file

    allocate(igp)                     !<-- Initialize hit_grid with defaults

    call igp%init(inputfile)          !<-- Properly initialize the hit_grid solver (see hit_grid.F90)
  
    call igp%start_io(.true.)                !<-- Start I/O by creating a header file (see io.F90)

    call igp%printDivergence()
    
    ! Fringe associations for non-periodic BCs in x
    call igp%fringe_x1%allocateTargetArray_Cells(utarget)                !<-- Allocate target array of appropriate size
    call igp%fringe_x1%allocateTargetArray_Cells(vtarget)                !<-- Allocate target array of appropriate size
    call igp%fringe_x1%allocateTargetArray_Edges(wtarget)                !<-- Allocate target array of appropriate size
    utarget = 1.d0                                                      !<-- Target u - velocity at inlet 
    vtarget = 0.d0                                                      !<-- Target v - velocity at inlet
    wtarget = 0.d0                                                      !<-- Target w - velocity at inlet    
    call igp%fringe_x1%associateFringeTargets(utarget, vtarget, wtarget) !<-- Link the target velocity array to igp 
    call igp%fringe_x2%associateFringeTargets(utarget, vtarget, wtarget) !<-- Link the target velocity array to igp 
  
    call tic() 
    do while (igp%tsim < igp%tstop) 
       
       call igp%timeAdvance()     !<-- Time stepping scheme + Pressure Proj. (see igrid.F90)
       call doTemporalStuff(igp)     !<-- Go to the temporal hook (see temporalHook.F90)
       
    end do 
 
    call igp%finalize_io()                  !<-- Close the header file (wrap up i/o)

    call igp%destroy()                !<-- Destroy the IGRID derived type 
   

    deallocate(igp)                   !<-- Deallocate all the memory associated with scalar defaults
    
    call MPI_Finalize(ierr)           !<-- Terminate MPI 

end program
