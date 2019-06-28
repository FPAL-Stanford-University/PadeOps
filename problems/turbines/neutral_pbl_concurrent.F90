! Template for PadeOps
#include "neutral_pbl_concurrent_files/initialize.F90"       
#include "neutral_pbl_concurrent_files/temporalHook.F90"  

program neutral_pbl_concurrent
    use mpi
    use kind_parameters,  only: clen, rkind
    use IncompressibleGrid, only: igrid
    use temporalhook, only: doTemporalStuff
    use timer, only: tic, toc
    use exits, only: message
    use constants, only: one, zero
    use budgets_time_avg_mod, only: budgets_time_avg  
    implicit none

    type(igrid), allocatable, target :: primary, precursor 
    character(len=clen) :: inputfile, primary_inputFile, precursor_inputFile
    integer :: ierr, ioUnit
    type(budgets_time_avg) :: budg_tavg
    real(rkind) :: dt1, dt2, dt

    namelist /concurrent/ primary_inputfile, precursor_inputfile

    call MPI_Init(ierr)                                                

    call GETARG(1,inputfile)                                            

    allocate(precursor)                                                       
    allocate(primary)                                                     
    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=concurrent)
    close(ioUnit)    

    call primary%init(primary_inputFile, .true.)                               
    call primary%start_io(.true.)                                          
    call primary%printDivergence()
   
    call precursor%init(precursor_inputFile, .false.)                                           
    precursor%Am_I_Primary = .false. 
    call precursor%start_io(.true.)                                           
   
    if (primary%usefringe) then
        call primary%fringe_x%associateFringeTargets(precursor%u, precursor%v, precursor%wC, precursor%T) 
    end if 

    call budg_tavg%init(primary_inputfile, primary)   !<-- Budget class initialization 

    call message("==========================================================")
    call message(0, "All memory allocated! Now running the simulation.")
    call tic() 
    do while (primary%tsim < primary%tstop) 
       dt1 = primary%get_dt(recompute=.true.)
       dt2 = precursor%get_dt(recompute=.true.)
       dt = min(dt1, dt2)
       
       call primary%timeAdvance(dt)

       call precursor%timeAdvance(dt)
       
       call budg_tavg%doBudgets()       

       call doTemporalStuff(primary,   1)                                        
       call doTemporalStuff(precursor,2)                                        
     
    end do 

    call budg_tavg%destroy()           !<-- release memory taken by the budget class 

    call precursor%finalize_io()          
    call primary%finalize_io()         

    call precursor%destroy()                
    call primary%destroy()               
   
    deallocate(precursor, primary)             
    
    call MPI_Finalize(ierr)        

end program
