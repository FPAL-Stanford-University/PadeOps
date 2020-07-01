! Template for PadeOps

#include "HIT_AD_interact_files/initialize.F90"       
#include "HIT_AD_interact_files/temporalHook.F90"  

program HIT_AD_interact
    use mpi
    use kind_parameters,  only: clen, rkind
    use IncompressibleGrid, only: igrid
    use temporalhook, only: doTemporalStuff
    use timer, only: tic, toc
    use exits, only: message, GracefulExit
    use constants, only: one, zero
    use HIT_AD_interact_parameters, only: simulationID, InflowProfileType, InflowProfileAmplit, InflowProfileThick
    use fof_mod, only: fof
    use budgets_time_avg_mod, only: budgets_time_avg  
    use budgets_vol_avg_mod, only: budgets_vol_avg  
    !use decomp_2d,                only: nrank
    implicit none

    type(igrid), allocatable, target :: hit, adsim
    character(len=clen) :: inputfile, HIT_InputFile, AD_InputFile, fof_dir, filoutdir
    integer :: ierr, ioUnit
    type(budgets_time_avg) :: budg_tavg
    type(budgets_vol_avg)  :: budg_vavg
    real(rkind), dimension(:,:,:), allocatable :: utarget0, vtarget0, wtarget0
    real(rkind), dimension(:,:,:), allocatable :: utarget1, vtarget1, wtarget1
    real(rkind) :: dt1, dt2, dt, InflowSpeed = 1.d0
    real(rkind) :: k_bandpass_left = 10.d0, k_bandpass_right = 64.d0, x_shift
    integer :: nxADSIM, nxHIT
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
    
    call adsim%fringe_x1%allocateTargetArray_Cells(utarget0)                
    call adsim%fringe_x1%allocateTargetArray_Cells(vtarget0)                
    call adsim%fringe_x1%allocateTargetArray_Edges(wtarget0)                
    select case(InflowProfileType)
    case(0)
        utarget0 = InflowSpeed
    case(1)
        !zinY => adsim%mesh(:,:,:,3)
        utarget0 = InflowSpeed*(one + InflowProfileAmplit*tanh((adsim%mesh(:,:,:,3)-adsim%zMid)/InflowProfileThick))
    end select
    vtarget0 = 0.d0                                                      
    wtarget0 = 0.d0                                                      
    call adsim%fringe_x1%associateFringeTargets(utarget0, vtarget0, wtarget0) 

    call adsim%fringe_x2%allocateTargetArray_Cells(utarget1)                
    call adsim%fringe_x2%allocateTargetArray_Cells(vtarget1)                
    call adsim%fringe_x2%allocateTargetArray_Edges(wtarget1)                
    utarget1 = 0.d0
    vtarget1 = 0.d0
    wtarget1 = 0.d0
    call adsim%fringe_x2%associateFringeTargets(utarget1, vtarget1, wtarget1) 

    nxHIT = hit%gpC%xsz(1)
    nxADSIM  = adsim%gpC%xsz(1)

    call hit%spectC%init_bandpass_filter(k_bandpass_left, k_bandpass_right, hit%cbuffzC(:,:,:,1), hit%cbuffyC(:,:,:,1))

    ! Set the true target field for AD simulation
    x_shift = adsim%tsim*InflowSpeed 
    call hit%spectC%bandpassFilter_and_phaseshift(hit%whatC , uTarget1(nxADSIM-nxHIT+1:nxADSIM,:,:), x_shift)
    call hit%interpolate_cellField_to_edgeField(uTarget1(nxADSIM-nxHIT+1:nxADSIM,:,:), wTarget1(nxADSIM-nxHIT+1:nxADSIM,:,:),0,0)
    call hit%spectC%bandpassFilter_and_phaseshift(hit%uhat  , uTarget1(nxADSIM-nxHIT+1:nxADSIM,:,:), x_shift)
    call hit%spectC%bandpassFilter_and_phaseshift(hit%vhat  , vTarget1(nxADSIM-nxHIT+1:nxADSIM,:,:), x_shift)
    
    ! Now scale rhw HIT field appropriately
    ! Note that the bandpass filtered velocity field has zero mean 
    uTarget1(nxADSIM-nxHIT+1:nxADSIM,:,:) = uTarget1(nxADSIM-nxHIT+1:nxADSIM,:,:)  + uTarget0(nxADSIM-nxHIT+1:nxADSIM,:,:)
    vTarget1(nxADSIM-nxHIT+1:nxADSIM,:,:) = vTarget1(nxADSIM-nxHIT+1:nxADSIM,:,:)
    wTarget1(nxADSIM-nxHIT+1:nxADSIM,:,:) = wTarget1(nxADSIM-nxHIT+1:nxADSIM,:,:)
   

    ! Graceful exit if apply filter is true and shear type is not 0
    if((InflowProfileType .ne.0) .and. applyfilters) then
        call GracefulExit("Scale separation on the fly is not supported for non-periodic in z. Set applyfilter to false", 111)
    endif

    ! Initialize the filters
    if (applyfilters) then
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

       call budg_tavg%doBudgets()       !<--- perform budget related operations -----Question::Where should this be placed ??
       call budg_vavg%doBudgets()       !<--- perform budget related operations -----Question::Where should this be placed ??

       x_shift = adsim%tsim*InflowSpeed
       call hit%spectC%bandpassFilter_and_phaseshift(hit%whatC , uTarget1(nxADSIM-nxHIT+1:nxADSIM,:,:), x_shift)
       call hit%interpolate_cellField_to_edgeField(uTarget1(nxADSIM-nxHIT+1:nxADSIM,:,:), wTarget1(nxADSIM-nxHIT+1:nxADSIM,:,:),0,0)
       call hit%spectC%bandpassFilter_and_phaseshift(hit%uhat  , uTarget1(nxADSIM-nxHIT+1:nxADSIM,:,:), x_shift)
       call hit%spectC%bandpassFilter_and_phaseshift(hit%vhat  , vTarget1(nxADSIM-nxHIT+1:nxADSIM,:,:), x_shift)
      
       ! Now scale rhw HIT field appropriately
       uTarget1(nxADSIM-nxHIT+1:nxADSIM,:,:) = uTarget1(nxADSIM-nxHIT+1:nxADSIM,:,:) + uTarget0(nxADSIM-nxHIT+1:nxADSIM,:,:)
       vTarget1(nxADSIM-nxHIT+1:nxADSIM,:,:) = vTarget1(nxADSIM-nxHIT+1:nxADSIM,:,:)
       wTarget1(nxADSIM-nxHIT+1:nxADSIM,:,:) = wTarget1(nxADSIM-nxHIT+1:nxADSIM,:,:)
     

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

    call budg_tavg%destroy()           !<-- release memory taken by the budget class 
    call budg_vavg%destroy()           !<-- release memory taken by the budget class 

    if (applyfilters) then
      do fid = 1,nfilters
         call filt(fid)%destroy()
      end do
      if (allocated(pid)) deallocate(pid)
      deallocate(filt)
    end if 

    call hit%finalize_io()          
    call adsim%finalize_io()         

    call hit%destroy()                
    call adsim%destroy()               
   
    deallocate(hit, adsim)             
    
    call MPI_Finalize(ierr)        

end program
