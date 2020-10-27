! Template for PadeOps

#include "concurrentSimulation_files/initialize.F90"       
#include "concurrentSimulation_files/temporalHook.F90"  

program concurrentSimulation
    use mpi
    use kind_parameters,  only: clen, rkind
    use IncompressibleGrid, only: igrid
    use temporalhook, only: doTemporalStuff
    use timer, only: tic, toc
    use exits, only: message, GracefulExit
    use budgets_time_avg_mod, only: budgets_time_avg  
    use budgets_xy_avg_mod, only: budgets_xy_avg  
    use decomp_2d, only: nrank
    use interpolatorMod, only: interpolator 
    use constants, only: zero, half

    implicit none

    type(igrid), allocatable, target :: igp, prec
    character(len=clen) :: inputfile, precInputFile, mainInputFile
    integer :: ierr
    integer :: ioUnit, prec_ixen_2
    type(interpolator) :: interpC, interpE
    type(budgets_time_avg) :: budg_tavg
    type(budgets_xy_avg) :: budg_xyavg
    logical :: useInterpolator = .false.
    real(rkind) :: srcxst, srcyst, srczst
    real(rkind), allocatable, dimension(:,:,:), target :: utarget, vtarget, wtarget, Ttarget
    namelist /concurrent/ precInputFile, mainInputFile, useInterpolator, srcxst, srcyst, srczst

    call MPI_Init(ierr)                                                 !<-- Begin MPI
    call GETARG(1,inputfile)                                            !<-- Get the location of the input file

    allocate(igp)                                                       !<-- Initialize hit_grid with defaults
    allocate(prec)                                                      !<-- Initialize hit_grid with defaults

    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=concurrent)
    close(ioUnit)    

    call prec%init(precInputFile, .true.)                                    !<-- Properly initialize the hit_grid solver (see hit_grid.F90)
    if(nrank==0) then
        print *, '---------- Done prec init---------  '
    endif
    call prec%start_io(.true.)                                           !<-- Start I/O by creating a header file (see io.F90)
    call prec%printDivergence()
 
    call igp%init(mainInputFile, .false.)                                            !<-- Properly initialize the hit_grid solver (see hit_grid.F90)
    if(nrank==0) then
        print *, '---------- Done main init---------  '
    endif
    call igp%start_io(.true.)                                           !<-- Start I/O by creating a header file (see io.F90)
    call igp%printDivergence()

    ! Fringe associations for non-periodic BCs in x
    if(useInterpolator) then
        ! Safeguards
        if( srcxst + igp%Lx > prec%Lx) then
            call GracefulExit("Precursor is not long enough. Check details.", 11)
        endif
        if( srcyst + igp%Ly > prec%Ly) then
            call GracefulExit("Precursor is not wide enough. Check details.", 11)
        endif
        if( srczst + igp%Lz > prec%Lz) then
            call GracefulExit("Precursor is not tall enough. Check details.", 11)
        endif
        if( abs(srczst) > zero) then
            call GracefulExit("srczst must be zero. Check details.", 11)
        endif
        if(.not. igp%useFringe) then
            call GracefulExit("If interpolator is being used, igp%useFringe must be true.", 11)
        endif
        !if( abs(srcyst+igp%Ly-prec%Ly) > 1.0d-12 ) then
        !  !if(.not. igp%fringe_x%Apply_y_fringe) then
        !  !    call GracefulExit("If part of the domain in y is being forced, igp%fringe_x%Apply_y_fringe must be true.", 11)
        !  !endif
        !endif

        ! Create the interpolator 
        write(*,*) "srcyst = ", srcyst
        write(*,*) "igp%yline(1) = ", igp%yline(1:3)
        call interpC%init(prec%gpC,igp%gpC,prec%xline,prec%yline,prec%zline, igp%xline+srcxst,igp%yline+srcyst,igp%zline +srczst, 'interpC_debug')
        call interpE%init(prec%gpE,igp%gpE,prec%xline,prec%yline,prec%zEline,igp%xline+srcxst,igp%yline+srcyst,igp%zEline+srczst,' interpE_debug')

        ! allocate memory
        allocate(utarget(igp%gpC%xsz(1),igp%gpC%xsz(2),igp%gpC%xsz(3)))
        allocate(vtarget(igp%gpC%xsz(1),igp%gpC%xsz(2),igp%gpC%xsz(3)))
        allocate(wtarget(igp%gpE%xsz(1),igp%gpE%xsz(2),igp%gpE%xsz(3)))
        if(igp%isStratified) then
            allocate(Ttarget(igp%gpE%xsz(1),igp%gpE%xsz(2),igp%gpE%xsz(3)))
        endif

        ! associate memory with pointers in igp%fringe_x
        if(igp%isStratified) then
            call igp%fringe_x%associateFringeTargets(utarget, vtarget, wtarget, Ttarget)
        else
            call igp%fringe_x%associateFringeTargets(utarget, vtarget, wtarget)
        endif

    else
        ! Link the target velocity array to igp 
        prec_ixen_2 = prec%fringe_x%get_ixen_2()
        call igp%fringe_x%associateFringeTargets(prec%u, prec%v, prec%w, xen_targ=prec_ixen_2)
    endif

    call budg_xyavg%init(precInputFile, prec)   !<-- Budget class initialization 

    call budg_tavg%init(mainInputFile, igp)   !<-- Budget class initialization 
    ! NOTE: Beyond this point, the target arrays can change in time, In case of
    ! concurrent simulations where inflow is turbulent, the utarget, vtarget and
    ! wtarget are simply set equal to the u, v, and w arrays of the concurrent
    ! simulation (at the end of each time step).  

    call tic() 
    do while (igp%tsim < igp%tstop) 
       call prec%timeAdvance()                                           !<- Time stepping scheme + Pressure Proj. (see igrid.F90)

       if(useInterpolator) then
           call interpC%LinInterp3D(prec%u, utarget)
           call interpC%LinInterp3D(prec%v, vtarget)
           call interpE%LinInterp3D(prec%w, wtarget)
           if(igp%isStratified) then
               call interpC%LinInterp3D(prec%T, Ttarget)
           endif
       endif

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
   
    deallocate(utarget, vtarget, wtarget)
    if(allocated(Ttarget)) deallocate(Ttarget)
    call interpE%destroy()
    call interpC%destroy()

    deallocate(prec)                                                     !<-- Deallocate all the memory associated with scalar defaults
    deallocate(igp)                                                      !<-- Deallocate all the memory associated with scalar defaults
    
    call MPI_Finalize(ierr)                                              !<-- Terminate MPI 

end program
