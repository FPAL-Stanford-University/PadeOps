! Template for PadeOps

#include "concurrentSimulationMult_files/initialize.F90"       
#include "concurrentSimulationMult_files/temporalHook.F90"  

program concurrentSimulationMult
    use mpi
    use kind_parameters,  only: clen, rkind
    use IncompressibleGrid, only: igrid
    use temporalhook, only: doTemporalStuff
    use timer, only: tic, toc
    use exits, only: message

    implicit none

    type(igrid), allocatable, target :: igp1, igp2, prec
    character(len=clen) :: inputfile, precInputFile, main1InputFile, main2InputFile
    integer :: ierr
    !real(rkind), dimension(:,:,:), allocatable :: utarget, vtarget, wtarget
    integer :: ioUnit
    namelist /concurrent/ precInputFile, main1InputFile, main2InputFile

    call MPI_Init(ierr)                                     !<-- Begin MPI
    call GETARG(1,inputfile)                                !<-- Get the location of the input file

    allocate(igp1)                                          !<-- Initialize hit_grid with defaults
    allocate(igp2)                                          !<-- Initialize hit_grid with defaults
    allocate(prec)                                          !<-- Initialize hit_grid with defaults

    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=concurrent)
    close(ioUnit)    

    call prec%init(precInputFile, .true.)                   !<-- Properly initialize the hit_grid solver (see hit_grid.F90)
    call prec%start_io(.true.)                              !<-- Start I/O by creating a header file (see io.F90)
    call prec%printDivergence()
 
    call igp1%init(main1InputFile, .false.)                 !<-- Properly initialize the hit_grid solver (see hit_grid.F90)
    call igp1%start_io(.true.)                              !<-- Start I/O by creating a header file (see io.F90)
    call igp1%printDivergence()

    call igp2%init(main2InputFile, .false.)                 !<-- Properly initialize the hit_grid solver (see hit_grid.F90)
    call igp2%start_io(.true.)                              !<-- Start I/O by creating a header file (see io.F90)
    call igp2%printDivergence()

    ! Fringe associations for non-periodic BCs in x
    call igp1%fringe_x%associateFringeTargets(prec%u, prec%v, prec%w) !<-- Link the target velocity array to igp 
    call igp2%fringe_x%associateFringeTargets(prec%u, prec%v, prec%w) !<-- Link the target velocity array to igp 

    ! NOTE: Beyond this point, the target arrays can change in time, In case of
    ! concurrent simulations where inflow is turbulent, the utarget, vtarget and
    ! wtarget are simply set equal to the u, v, and w arrays of the concurrent
    ! simulation (at the end of each time step).  

    call tic() 
    do while (igp2%tsim < igp2%tstop) 
       
       call prec%timeAdvance()                                   !<- Time stepping scheme + Pressure Proj. (see igrid.F90)
       call igp1%timeAdvance(prec%get_dt())                      !<-- Time stepping scheme + Pressure Proj. (see igrid.F90)
       call igp2%timeAdvance(prec%get_dt())                      !<-- Time stepping scheme + Pressure Proj. (see igrid.F90)
       call doTemporalStuff(prec, igp1, igp2)                    !<-- Go to the temporal hook (see temporalHook.F90)
       
    end do 
 
    call prec%finalize_io()                                      !<-- Close the header file (wrap up i/o)
    call igp1%finalize_io()                                      !<-- Close the header file (wrap up i/o)
    call igp2%finalize_io()                                      !<-- Close the header file (wrap up i/o)

    call prec%destroy()                                          !<-- Destroy the IGRID derived type 
    call igp1%destroy()                                          !<-- Destroy the IGRID derived type 
    call igp2%destroy()                                          !<-- Destroy the IGRID derived type 
   

    deallocate(prec)                                             !<-- Deallocate all the memory associated with scalar defaults
    deallocate(igp1)                                             !<-- Deallocate all the memory associated with scalar defaults
    deallocate(igp2)                                             !<-- Deallocate all the memory associated with scalar defaults
    
    call MPI_Finalize(ierr)                                      !<-- Terminate MPI 

end program
