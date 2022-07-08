! Template for PadeOps

#include "HIT_AD_interact_files/initialize.F90"       
#include "HIT_AD_interact_files/temporalHook.F90"  

program HIT_AD_interact
    use mpi
    use kind_parameters,            only: clen, rkind
    use IncompressibleGrid,         only: igrid
    use temporalhook,               only: doTemporalStuff
    use timer,                      only: tic, toc
    use exits,                      only: message, GracefulExit
    use fortran_assert,             only: assert
    use constants,                  only: one, zero
    use HIT_AD_interact_parameters, only: simulationID, InflowProfileType, &
                                          InflowProfileAmplit, InflowProfileThick, &
                                          streamWiseCoord, nxHIT, nyHIT, nzHIT, &
                                          nxADSIM, nyADSIM, nzADSIM, copyHITfieldsToADSIM
    use fof_mod,                    only: fof
    use budgets_time_avg_mod,       only: budgets_time_avg  
    use budgets_vol_avg_mod,        only: budgets_vol_avg  
    !use decomp_2d,                only: nrank
    implicit none

    type(igrid), allocatable, target :: hit, adsim
    character(len=clen) :: inputfile, HIT_InputFile, AD_InputFile, fof_dir, filoutdir
    integer :: ierr, ioUnit
    type(budgets_time_avg) :: budg_tavg
    type(budgets_vol_avg)  :: budg_vavg
    real(rkind), dimension(:,:,:), allocatable :: utarget0, vtarget0, wtarget0
    real(rkind), dimension(:,:,:), allocatable :: utarget1, vtarget1, wtarget1
    real(rkind), dimension(:,:,:), allocatable :: uhitFilt, vhitFilt, whitFilt
    
    real(rkind) :: dt1, dt2, dt, InflowSpeed = 1.d0
    real(rkind) :: k_bandpass_left = 10.d0, k_bandpass_right = 64.d0, x_shift
    type(fof), dimension(:), allocatable :: filt
    integer, dimension(:), allocatable :: pid
    integer :: fid, nfilters = 2, tid_FIL_FullField = 75, tid_FIL_Planes = 4
    logical :: applyFilters = .false. 
    logical, parameter :: synchronize_RK_substeps = .true.

    namelist /concurrent/ HIT_InputFile, AD_InputFile, InflowSpeed, k_bandpass_left, k_bandpass_right
    namelist /FILTER_INFO/ applyfilters, nfilters, fof_dir, tid_FIL_FullField, tid_FIL_Planes, filoutdir

    call MPI_Init(ierr)                                                

    call GETARG(1,inputfile)                                            

    allocate(hit)                                                       
    allocate(adsim)                                                     
    
    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=concurrent)
    read(unit=ioUnit, NML=FILTER_INFO)
    close(ioUnit)    

    simulationID = 1
    call adsim%init(AD_InputFile, .true.)                               
    call adsim%start_io(.true.)                                          
    call adsim%printDivergence()
   
    call mpi_barrier(mpi_comm_world, ierr)

    simulationID = 2
    call hit%init(HIT_InputFile, .false.)                                           
    hit%Am_I_Primary = .false. 
    call hit%start_io(.true.)                                           
    call hit%printDivergence()

    allocate(uhitFilt(hit%gpC%xsz(1),hit%gpC%xsz(2),hit%gpC%xsz(3)))
    allocate(vhitFilt(hit%gpC%xsz(1),hit%gpC%xsz(2),hit%gpC%xsz(3)))
    allocate(whitFilt(hit%gpE%xsz(1),hit%gpE%xsz(2),hit%gpE%xsz(3)))
    
    call adsim%fringe_x1%allocateTargetArray_Cells(utarget0)                
    call adsim%fringe_x1%allocateTargetArray_Cells(vtarget0)                
    call adsim%fringe_x1%allocateTargetArray_Edges(wtarget0)                
    select case(InflowProfileType)
    case(0)
      select case (streamWiseCoord)
        case ('x')
          utarget0 = InflowSpeed
          vtarget0 = 0.d0                                                      
          wtarget0 = 0.d0                                                      
        case ('y')
          utarget0 = 0.d0                                                      
          vtarget0 = InflowSpeed
          wtarget0 = 0.d0                                                      
        case ('z')
          utarget0 = 0.d0                                                      
          vtarget0 = 0.d0                                                      
          wtarget0 = InflowSpeed
      end select
    case(1)
      select case (streamWiseCoord)
        case ('x')
          !zinY => adsim%mesh(:,:,:,3)
          utarget0 = InflowSpeed*(one + InflowProfileAmplit*tanh((adsim%mesh(:,:,:,3)-adsim%zMid)/InflowProfileThick))
        case ('y')
          call gracefulExit('The code only supports x-direction streamwise'//&
            'coordinate at this time',123)
        case ('z')
          call gracefulExit('The code only supports x-direction streamwise'//&
            'coordinate at this time',123)
      end select
    end select
    call adsim%fringe_x1%associateFringeTargets(utarget0, vtarget0, wtarget0) 

    call adsim%fringe_x2%allocateTargetArray_Cells(utarget1)                
    call adsim%fringe_x2%allocateTargetArray_Cells(vtarget1)                
    call adsim%fringe_x2%allocateTargetArray_Edges(wtarget1)                
    utarget1 = 0.d0
    vtarget1 = 0.d0
    wtarget1 = 0.d0
    call adsim%fringe_x2%associateFringeTargets(utarget1, vtarget1, wtarget1) 

    nxHIT = hit%gpC%xsz(1)
    nyHIT = hit%gpC%ysz(2)
    nzHIT = hit%gpC%zsz(3)
    nxADSIM  = adsim%gpC%xsz(1)
    nyADSIM  = adsim%gpC%ysz(2)
    nzADSIM  = adsim%gpC%zsz(3)

    call hit%spectC%init_bandpass_filter(k_bandpass_left, k_bandpass_right, hit%cbuffzC(:,:,:,1), hit%cbuffyC(:,:,:,1))

  

    ! Set the true target field for AD simulation
    x_shift = adsim%tsim*InflowSpeed 
    call hit%spectC%bandpassFilter_and_phaseshift(hit%whatC , uhitFilt, x_shift, streamWiseCoord)
    call hit%interpolate_cellField_to_edgeField(uhitFilt,whitFilt,0,0)
    call hit%spectC%bandpassFilter_and_phaseshift(hit%uhat  , uhitFilt, x_shift, streamWiseCoord)
    call hit%spectC%bandpassFilter_and_phaseshift(hit%vhat  , vhitFilt, x_shift, streamWiseCoord)

    ! Go from hitFilt to uTargets  
    call copyHITfieldsToADSIM(uhitFilt,vhitFilt,whitFilt,utarget1,vtarget1,wtarget1,hit,adsim,streamWiseCoord)

    ! Now scale rhw HIT field appropriately
    ! Note that the bandpass filtered velocity field has zero mean 
    uTarget1 = uTarget1  + uTarget0
    vTarget1 = vTarget1  + vTarget0
    wTarget1 = wTarget1  + wTarget0
   

    ! Graceful exit if apply filter is true and shear type is not 0
    if((InflowProfileType .ne.0) .and. applyfilters) then
        call GracefulExit("Scale separation on the fly is not supported for non-periodic in z. Set applyfilter to false", 111)
    endif

    ! Initialize the filters
    if (applyfilters) then
      call assert(streamWiseCoord == 'x','Filters only supported for '//&
        'streamwise coordinate in x')
      allocate(filt(nfilters))
      do fid = 1,nfilters
         call filt(fid)%init(adsim%runID, fof_dir, filoutdir, fid, adsim%spectC, adsim%cbuffyC(:,:,:,1), adsim%cbuffzC(:,:,:,1), .true., adsim%gpC)
      end do 
      if (allocated(adsim%yplanes)) then 
         allocate(pid(size(adsim%yplanes)))
         pid = adsim%yplanes 
      end if 
    end if 

    call budg_tavg%init(AD_Inputfile, adsim)   !<-- Budget class initialization 
    call budg_vavg%init(HIT_Inputfile, hit)    !<-- Budget class initialization 

    call message("==========================================================")
    call message(0, "All memory allocated! Now running the simulation.")
    call tic() 
    do while (adsim%tsim < adsim%tstop) 
       dt1 = adsim%get_dt(recompute=.true.)
       dt2 = hit%get_dt(recompute=.true.)
       dt = min(dt1, dt2)
      
       if (synchronize_RK_substeps) then
           adsim%dt = dt
           hit%dt = dt
           ! Stage 1
           call adsim%advance_SSP_RK45_Stage_1()
           call hit%advance_SSP_RK45_Stage_1()
           ! Stage 2
           call adsim%advance_SSP_RK45_Stage_2()
           call hit%advance_SSP_RK45_Stage_2()
           ! Stage 3
           call adsim%advance_SSP_RK45_Stage_3()
           call hit%advance_SSP_RK45_Stage_3()
           ! Stage 4
           call adsim%advance_SSP_RK45_Stage_4()
           call hit%advance_SSP_RK45_Stage_4()
           ! Stage 5
           call adsim%advance_SSP_RK45_Stage_5()
           call hit%advance_SSP_RK45_Stage_5()
           ! Call wrap up 
           call adsim%wrapup_timestep()
           call hit%wrapup_timestep() 

       else
          call adsim%timeAdvance(dt)
          call hit%timeAdvance(dt)
       end if  

       call budg_tavg%doBudgets()       !<--- perform budget related operations
       call budg_vavg%doBudgets()       !<--- perform budget related operations

       x_shift = adsim%tsim*InflowSpeed
       call hit%spectC%bandpassFilter_and_phaseshift(hit%whatC , &
        uhitFilt, x_shift, streamWiseCoord)
       call hit%interpolate_cellField_to_edgeField(uhitFilt, &
         whitFilt,0,0)
       call hit%spectC%bandpassFilter_and_phaseshift(hit%uhat  , &
         uhitFilt, x_shift, streamWiseCoord)
       call hit%spectC%bandpassFilter_and_phaseshift(hit%vhat  , &
         vhitFilt, x_shift, streamWiseCoord)


       ! Now copy hitFilt into Target1 
       call copyHITfieldsToADSIM(uhitFilt,vhitFilt,whitFilt,utarget1,vtarget1,wtarget1,hit,adsim,streamWiseCoord)

       ! Now scale rhw HIT field appropriately
       uTarget1 = uTarget1 + uTarget0 
       vTarget1 = vTarget1 + vTarget0 
       wTarget1 = wTarget1 + wTarget0 

     

       call doTemporalStuff(adsim, 1)                                        
       call doTemporalStuff(hit  , 2)                                        
     

       if (applyfilters) then
         if ((mod(adsim%step,tid_FIL_FullField)==0) .or. (mod(adsim%step,tid_FIL_Planes)==0))  then
            do fid = 1,nfilters
               call filt(fid)%filter_Complex2Real(adsim%uhat, adsim%rbuffxC(:,:,:,1))
               if (mod(adsim%step,tid_FIL_FullField) == 0) call filt(fid)%dumpFullField( adsim%rbuffxC(:,:,:,1), "uVel",adsim%step)
               if (mod(adsim%step,tid_FIL_Planes) == 0) call filt(fid)%dumpYplanes( adsim%rbuffxC(:,:,:,1), "uVel",pid, adsim%step)
               
               call filt(fid)%filter_Complex2Real(adsim%vhat, adsim%rbuffxC(:,:,:,1))
               if (mod(adsim%step,tid_FIL_FullField) == 0) call filt(fid)%dumpFullField( adsim%rbuffxC(:,:,:,1), "vVel",adsim%step)
               if (mod(adsim%step,tid_FIL_Planes) == 0) call filt(fid)%dumpYplanes( adsim%rbuffxC(:,:,:,1), "vVel",pid, adsim%step)

               call filt(fid)%filter_Complex2Real(adsim%whatC, adsim%rbuffxC(:,:,:,1))
               if (mod(adsim%step,tid_FIL_FullField) == 0) call filt(fid)%dumpFullField( adsim%rbuffxC(:,:,:,1), "wVel",adsim%step)
               if (mod(adsim%step,tid_FIL_Planes) == 0) call filt(fid)%dumpYplanes( adsim%rbuffxC(:,:,:,1), "wVel",pid, adsim%step)

               if (allocated(adsim%Pressure)) then
                  call filt(fid)%filter_Real2Real(adsim%pressure, adsim%rbuffxC(:,:,:,1))
                  if (mod(adsim%step,tid_FIL_FullField) == 0) call filt(fid)%dumpFullField( adsim%rbuffxC(:,:,:,1), "pres",adsim%step)
                  if (mod(adsim%step,tid_FIL_Planes) == 0) call filt(fid)%dumpYplanes( adsim%rbuffxC(:,:,:,1), "pres",pid, adsim%step)
               end if
            end do
            call message(0,"Just dumped filtered data")
         end if
       end if

    end do

    call budg_tavg%doBudgets(.true.)   !<--- force dump if budget calculation had started
    call budg_vavg%doBudgets(.true.)   !<--- force dump if budget calculation had started

    call budg_tavg%destroy()           !<-- release memory taken by the budget class
    call budg_vavg%destroy()           !<-- release memory taken by the budget class

    if (applyfilters) then
      do fid = 1,nfilters
         call filt(fid)%destroy()
      end do
      if (allocated(pid)) deallocate(pid)
      deallocate(filt)
    end if

    deallocate(uhitFilt, vhitFilt, whitFilt)
    deallocate(utarget0, vtarget0, wtarget0)
    deallocate(utarget1, vtarget1, wtarget1)
    
    call hit%finalize_io()
    call adsim%finalize_io()

    call hit%destroy()
    call adsim%destroy()
   
    deallocate(hit, adsim)
    
    call MPI_Finalize(ierr)

end program
