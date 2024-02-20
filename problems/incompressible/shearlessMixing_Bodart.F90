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
    use stats_xy_mod,       only: stats_xy
    use timer,              only: tic, toc 

    implicit none

    type(igrid), allocatable, target :: SM
    type(stats_xy) :: stats
    character(len=clen) :: inputfile
    integer :: ierr
    
    call MPI_Init(ierr)

    call GETARG(1,inputfile)                                            

    allocate(SM)                                                     
    
    call SM%init(inputfile, .true.)                              
    call SM%start_io(.true.)                                          
    call SM%printDivergence()
    call stats%init(inputfile,SM)

    ! Get t=0 stats
    !call stats%compute_stats()

    call message("==========================================================")
    call message(0, "All memory allocated! Now running the simulation.")
    do while (SM%tsim < SM%tstop .and. SM%step < SM%nsteps)
       call tic() 
       call SM%timeAdvance()
       call stats%compute_stats()
       call doTemporalStuff(SM) 
       call toc()
    end do

    call message("==========================================================")
    call message(0,"Finalizing simulation")
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    call stats%destroy()
    call SM%finalize_io()
    call SM%destroy()
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    call message(0,"Simulation finalized")
   
    !deallocate(SM)
    
    call MPI_Finalize(ierr)

end program
