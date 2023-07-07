! Template for PadeOps

#include "shearLayerEntrainment_files/initialize.F90"       
#include "shearLayerEntrainment_files/temporalHook.F90"  

program shearLayerEntrainment
    use mpi
    use kind_parameters,    only: clen, rkind
    use IncompressibleGrid, only: igrid
    use temporalhook,       only: doTemporalStuff
    use exits,              only: message
    use budgets_xy_avg_mod, only: budgets_xy_avg
    use shearLayerEntrainment_parameters, only: initializeForceLayer 

    implicit none

    type(igrid), allocatable, target :: sim
    type(budgets_xy_avg) :: budg
    character(len=clen) :: inputfile
    integer :: ierr
    
    call MPI_Init(ierr)                                                

    call GETARG(1,inputfile)                                            

    allocate(sim)                                                     
    
    call sim%init(inputfile, .true.)
    call initializeForceLayer(sim) 
    call sim%start_io(.true.)                                          
    call sim%printDivergence()
    call budg%init(inputfile,sim)

    call message("==========================================================")
    call message(0, "All memory allocated! Now running the simulation.")
    do while (sim%tsim < sim%tstop .and. sim%step < sim%nsteps) 
       call sim%timeAdvance()
       call budg%doBudgets()
       call doTemporalStuff(sim) 
    end do

    call sim%finalize_io()
    call budg%destroy()
    call sim%destroy()
   
    !deallocate(sim)
    
    call MPI_Finalize(ierr)

end program
