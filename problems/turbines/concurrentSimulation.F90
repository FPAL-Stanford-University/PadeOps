! Template for PadeOps

#include "concurrentSimulation_files/initialize.F90"       
#include "concurrentSimulation_files/temporalHook.F90"  

program concurrentSimulation
    use mpi
    use kind_parameters,  only: clen, rkind
    use IncompressibleGrid, only: igrid
    use temporalhook, only: doTemporalStuff
    use timer, only: tic, toc
    use exits, only: message
    use budgets_time_avg_mod, only: budgets_time_avg  
    use budgets_xy_avg_mod, only: budgets_xy_avg  

    implicit none

    type(igrid), allocatable, target :: igp, prec
    character(len=clen) :: inputfile, precInputFile, mainInputFile
    integer :: ierr
    !real(rkind), dimension(:,:,:), allocatable :: utarget, vtarget, wtarget
    integer :: ioUnit
    type(budgets_time_avg) :: budg_tavg
    type(budgets_xy_avg) :: budg_xyavg
    namelist /concurrent/ precInputFile, mainInputFile

    call MPI_Init(ierr)                                                 !<-- Begin MPI
    call GETARG(1,inputfile)                                            !<-- Get the location of the input file

    allocate(igp)                                                       !<-- Initialize hit_grid with defaults
    allocate(prec)                                                      !<-- Initialize hit_grid with defaults

    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=concurrent)
    close(ioUnit)    

    call prec%init(precInputFile, .true.)                                    !<-- Properly initialize the hit_grid solver (see hit_grid.F90)
    call prec%start_io(.true.)                                           !<-- Start I/O by creating a header file (see io.F90)
    call prec%printDivergence()
 
    call igp%init(mainInputFile, .false.)                                            !<-- Properly initialize the hit_grid solver (see hit_grid.F90)
    call igp%start_io(.true.)                                           !<-- Start I/O by creating a header file (see io.F90)
    call igp%printDivergence()

    ! Fringe associations for non-periodic BCs in x
    call igp%fringe_x%associateFringeTargets(prec%u, prec%v, prec%w) !<-- Link the target velocity array to igp 

    call budg_xyavg%init(precInputFile, prec)   !<-- Budget class initialization 

    call budg_tavg%init(mainInputFile, igp)   !<-- Budget class initialization 
    ! NOTE: Beyond this point, the target arrays can change in time, In case of
    ! concurrent simulations where inflow is turbulent, the utarget, vtarget and
    ! wtarget are simply set equal to the u, v, and w arrays of the concurrent
    ! simulation (at the end of each time step).  

    call tic() 
    do while (igp%tsim < igp%tstop) 
       call prec%timeAdvance()                                           !<- Time stepping scheme + Pressure Proj. (see igrid.F90)
       call igp%timeAdvance(prec%get_dt())                               !<-- Time stepping scheme + Pressure Proj. (see igrid.F90)
       call budg_xyavg%doBudgets()       
       call budg_tavg%doBudgets()       
       call doTemporalStuff(prec, igp)                                   !<-- Go to the temporal hook (see temporalHook.F90)
       
    end do 

    call budg_xyavg%doBudgets(.true.)
    call budg_tavg%doBudgets(.true.)

    call budg_xyavg%destroy()
    call budg_tavg%destroy()
 
    call prec%finalize_io()                                              !<-- Close the header file (wrap up i/o)
    call igp%finalize_io()                                               !<-- Close the header file (wrap up i/o)

    call prec%destroy()                                                  !<-- Destroy the IGRID derived type 
    call igp%destroy()                                                   !<-- Destroy the IGRID derived type 
   

    deallocate(prec)                                                     !<-- Deallocate all the memory associated with scalar defaults
    deallocate(igp)                                                      !<-- Deallocate all the memory associated with scalar defaults
    
    call MPI_Finalize(ierr)                                              !<-- Terminate MPI 

end program
