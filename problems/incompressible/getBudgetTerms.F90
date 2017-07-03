! Template for PadeOps

#include "getBudgetTerms_files/initialize.F90"       
#include "getBudgetTerms_files/temporalHook.F90"  

program getBudgetTerms
    use mpi
    use kind_parameters,  only: clen
    use IncompressibleGrid, only: igrid
    !use IncompressibleGridWallM, only: igridWallM
    use temporalhook, only: doTemporalStuff
    use timer, only: tic, toc
    use exits, only: message
    use reductions, only: p_maxval
    use getBudgetTerms_parameters

    implicit none

    type(igrid), allocatable, target :: igp
    character(len=clen) :: inputfile
    integer :: ierr, tid

    call MPI_Init(ierr)               !<-- Begin MPI

    call GETARG(1,inputfile)          !<-- Get the location of the input file

    allocate(igp)                     !<-- Initialize hit_grid with defaults

    call igp%init(inputfile)          !<-- Properly initialize the hit_grid solver (see hit_grid.F90)
  
    !call igp%start_io(.true.)                !<-- Start I/O by creating a header file (see io.F90)

    !call igp%printDivergence()
  
    call tic() 
    !do while ((igp%tsim < igp%tstop) .and. (igp%step < igp%nsteps))
      
    do tid = tStatStart, tStatEnd, dtStats

       call message(0,"Time Step = ", tid)
       ! Ensure igp%runID is set correctly in setup
       igp%step = tid
       call igp%readFullField(igp%u, 'uVel')!; call message(1,"umax: ", p_maxval(igp%u))
       call igp%readFullField(igp%v, 'vVel')!; call message(1,"vmax: ", p_maxval(igp%v))
       call igp%readFullField(igp%w, 'wVel')!; call message(1,"wmax: ", p_maxval(igp%w))
       if (igp%isStratified .or. igp%initspinup) call igp%readFullField(igp%T,'potT')
       if (igp%fastCalcPressure) call igp%readFullField(igp%pressure,'prss')
     
        call igp%spectC%fft(igp%u,igp%uhat)   
        call igp%spectC%fft(igp%v,igp%vhat)   
        call igp%spectE%fft(igp%w,igp%what)   
        if (igp%isStratified .or. igp%initspinup) call igp%spectC%fft(igp%T,igp%That)   

       call igp%project_and_prep(igp%fastCalcPressure) 
       call igp%populate_rhs() 
       call igp%compute_stats3D()
 
       if (mod(tid,dtStatsDump)==0) then
         call message(0,"Dumping stats3D", tid)
         call igp%dump_stats3D()
       endif
    end do 
 
    call igp%finalize_io()                  !<-- Close the header file (wrap up i/o)

    call igp%destroy()                !<-- Destroy the IGRID derived type 
   

    deallocate(igp)                   !<-- Deallocate all the memory associated with scalar defaults
    
    call MPI_Finalize(ierr)           !<-- Terminate MPI 

end program
