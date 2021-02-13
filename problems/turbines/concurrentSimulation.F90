! Template for PadeOps

#include "concurrentSimulation_files/initialize.F90"       
#include "concurrentSimulation_files/temporalHook.F90"  

program concurrentSimulation
    use mpi
    use kind_parameters,  only: clen, rkind
    use IncompressibleGrid, only: igrid
    use temporalhook, only: doTemporalStuff, doTemporalStuff_prec, doTemporalStuff_gp
    use timer, only: tic, toc
    use exits, only: message, GracefulExit
    use budgets_time_avg_mod, only: budgets_time_avg  
    use budgets_xy_avg_mod, only: budgets_xy_avg  
    use decomp_2d, only: nrank
    use interpolatorMod, only: interpolator 
    use constants, only: zero, half, eps

    implicit none

    type(igrid), allocatable, target :: igp, prec
    character(len=clen) :: inputfile, precInputFile, mainInputFile
    integer :: ierr
    integer :: ioUnit, prec_ixen_2
    type(interpolator) :: interpC, interpE
    type(budgets_time_avg) :: budg_tavg
    type(budgets_xy_avg) :: budg_xyavg1, budg_xyavg2
    logical :: useInterpolator = .false., restrict_igp_dt = .false.
    real(rkind) :: srcxst, srcyst, srczst, tsub, dtfull, dtpart
    real(rkind), allocatable, dimension(:,:,:), target :: utarget,  vtarget,  wtarget,  Ttarget
    real(rkind), allocatable, dimension(:,:,:), target :: utarget0, vtarget0, wtarget0, Ttarget0
    real(rkind), allocatable, dimension(:,:,:), target :: utarget1, vtarget1, wtarget1, Ttarget1
    namelist /concurrent/ precInputFile, mainInputFile, useInterpolator, srcxst, srcyst, srczst, restrict_igp_dt

    call MPI_Init(ierr)                                                 !<-- Begin MPI
    call GETARG(1,inputfile)                                            !<-- Get the location of the input file

    allocate(igp)                                                       !<-- Initialize hit_grid with defaults
    allocate(prec)                                                      !<-- Initialize hit_grid with defaults

    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=concurrent)
    close(ioUnit)    
    call message(1, 'Starting precinit')

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

        !! Create the interpolator 
        !write(*,*) "srcyst = ", srcyst
        !write(*,*) "igp%yline(1) = ", igp%yline(1:3)
        call interpC%init(prec%gpC,igp%gpC,prec%xline,prec%yline,prec%zline, igp%xline+srcxst,igp%yline+srcyst,igp%zline +srczst, 'interpC_debug', igp%sgsmodel%get_z0())
        call interpE%init(prec%gpE,igp%gpE,prec%xline,prec%yline,prec%zEline,igp%xline+srcxst,igp%yline+srcyst,igp%zEline+srczst,' interpE_debug', igp%sgsmodel%get_z0())

        ! allocate memory
        allocate(utarget(igp%gpC%xsz(1),igp%gpC%xsz(2),igp%gpC%xsz(3)))
        allocate(vtarget(igp%gpC%xsz(1),igp%gpC%xsz(2),igp%gpC%xsz(3)))
        allocate(wtarget(igp%gpE%xsz(1),igp%gpE%xsz(2),igp%gpE%xsz(3)))
        allocate(utarget0(igp%gpC%xsz(1),igp%gpC%xsz(2),igp%gpC%xsz(3)))
        allocate(vtarget0(igp%gpC%xsz(1),igp%gpC%xsz(2),igp%gpC%xsz(3)))
        allocate(wtarget0(igp%gpE%xsz(1),igp%gpE%xsz(2),igp%gpE%xsz(3)))
        allocate(utarget1(igp%gpC%xsz(1),igp%gpC%xsz(2),igp%gpC%xsz(3)))
        allocate(vtarget1(igp%gpC%xsz(1),igp%gpC%xsz(2),igp%gpC%xsz(3)))
        allocate(wtarget1(igp%gpE%xsz(1),igp%gpE%xsz(2),igp%gpE%xsz(3)))
        if(igp%isStratified) then
            allocate(Ttarget(igp%gpE%xsz(1),igp%gpE%xsz(2),igp%gpE%xsz(3)))
            allocate(Ttarget0(igp%gpE%xsz(1),igp%gpE%xsz(2),igp%gpE%xsz(3)))
            allocate(Ttarget1(igp%gpE%xsz(1),igp%gpE%xsz(2),igp%gpE%xsz(3)))
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

    call budg_xyavg1%init(precInputFile, prec)   !<-- Budget class initialization 

    call budg_xyavg2%init(mainInputFile, igp)   !<-- Budget class initialization 
    call budg_tavg%init(mainInputFile, igp)   !<-- Budget class initialization 
    ! NOTE: Beyond this point, the target arrays can change in time, In case of
    ! concurrent simulations where inflow is turbulent, the utarget, vtarget and
    ! wtarget are simply set equal to the u, v, and w arrays of the concurrent
    ! simulation (at the end of each time step).  

    if(useInterpolator) then
        call interpC%LinInterp3D(prec%u, utarget1, .true.)
        call interpC%LinInterp3D(prec%v, vtarget1, .true.)
        !call interpC%loglaw_correction_uv(prec%u, prec%v, utarget1, vtarget1)
        call interpE%LinInterp3D(prec%w, wtarget1, .false.)
        if(igp%isStratified) then
            call interpC%LinInterp3D(prec%T, Ttarget1, .false.)
            !call interpC%loglaw_correction_T(prec%T, Ttarget1)
        endif
    endif

    call tic() 
    do while (igp%tsim < igp%tstop) 
       call prec%timeAdvance()                                           !<- Time stepping scheme + Pressure Proj. (see igrid.F90)
       call budg_xyavg1%doBudgets()       

       if(useInterpolator) then
           call doTemporalStuff_prec(prec)                                   !<-- Go to the temporal hook (see temporalHook.F90)
      
           utarget0 = utarget1; vtarget0 = vtarget1; wtarget0 = wtarget1; 
           if(igp%isStratified) Ttarget0 = Ttarget1
 
           call interpC%LinInterp3D(prec%u, utarget1, .true.)
           call interpC%LinInterp3D(prec%v, vtarget1, .true.)
           call interpE%LinInterp3D(prec%w, wtarget1, .false.)
           if(igp%isStratified) then
               call interpC%LinInterp3D(prec%T, Ttarget1, .false.)
           endif

           dtfull = prec%tsim-igp%tsim
           tsub = 0.0_rkind
           do while (igp%tsim < prec%tsim)

             dtpart = igp%get_dt(.true.)
             if(restrict_igp_dt) then
                 dtpart = min(dtpart, prec%tsim-igp%tsim+eps)         ! Ensure that igp does not outpace prec
             endif
             tsub = tsub + dtpart
             utarget = utarget0 + tsub/dtfull * (utarget1 - utarget0)
             vtarget = vtarget0 + tsub/dtfull * (vtarget1 - vtarget0)
             wtarget = wtarget0 + tsub/dtfull * (wtarget1 - wtarget0)
             if(igp%isStratified) then
               Ttarget = Ttarget0 + tsub/dtfull * (Ttarget1 - Ttarget0)
             endif

             call igp%timeAdvance(dtpart)                    !<-- Time stepping scheme + Pressure Proj. (see igrid.F90)
             call budg_xyavg2%doBudgets()       
             call budg_tavg%doBudgets()       
             call doTemporalStuff_gp(igp)                    !<-- Go to the temporal hook (see temporalHook.F90)
       
           end do 
       else
           call igp%timeAdvance(prec%get_dt())              !<-- Time stepping scheme + Pressure Proj. (see igrid.F90)
           call budg_xyavg2%doBudgets()       
           call budg_tavg%doBudgets()       
           call doTemporalStuff(prec, igp)                    !<-- Go to the temporal hook (see temporalHook.F90)
       endif

       !call igp%do_debug()

    end do 

    call budg_xyavg1%doBudgets(.true.)
    call budg_xyavg2%doBudgets(.true.)
    call budg_tavg%doBudgets(.true.)

    call budg_xyavg1%destroy()
    call budg_xyavg2%destroy()
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
