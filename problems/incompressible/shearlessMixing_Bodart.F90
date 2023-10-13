! Template for PadeOps

#include "shearlessMixing_Bodart_files/initialize.F90"       
#include "shearlessMixing_Bodart_files/temporalHook.F90"  

program shearlessMixing
    use mpi
    use kind_parameters,    only: clen, rkind
    use IncompressibleGrid, only: igrid
    use temporalhook,       only: doTemporalStuff
    use exits,              only: message
    use budgets_xy_avg_mod, only: budgets_xy_avg
    use timer,              only: tic, toc 

    implicit none

    type(igrid), allocatable, target :: SM
    type(budgets_xy_avg) :: budg
    character(len=clen) :: inputfile
    integer :: ierr
    
    call MPI_Init(ierr)

    call GETARG(1,inputfile)                                            

    allocate(SM)                                                     
    
    call SM%init(inputfile, .true.)                              
    call SM%start_io(.true.)                                          
    call SM%printDivergence()
    call budg%init(inputfile,SM)

    call message("==========================================================")
    call message(0, "All memory allocated! Now running the simulation.")
    do while (SM%tsim < SM%tstop .and. SM%step < SM%nsteps)
       call tic() 
       call SM%timeAdvance()
       call budg%doBudgets()
       call doTemporalStuff(SM) 
       call toc()
    end do

    call message("==========================================================")
    call message(0,"Finalizing simulation")
    call SM%finalize_io()
    call budg%destroy()
    call SM%destroy()
    call message(0,"Simulation finalized")
   
    !deallocate(SM)
    
    call MPI_Finalize(ierr)

end program
