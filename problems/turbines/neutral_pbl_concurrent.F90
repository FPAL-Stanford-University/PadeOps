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
    use HIT_AD_interact_parameters, only: simulationID
    use fof_mod, only: fof
    use budgets_time_avg_mod, only: budgets_time_avg  
    !use decomp_2d,                only: nrank
    implicit none

    type(igrid), allocatable, target :: hit, adsim
    character(len=clen) :: inputfile, HIT_InputFile, AD_InputFile, fof_dir, filoutdir
    integer :: ierr, ioUnit
    type(budgets_time_avg) :: budg_tavg
    real(rkind), dimension(:,:,:), allocatable :: utarget0, vtarget0, wtarget0
    real(rkind), dimension(:,:,:), allocatable :: utarget1, vtarget1, wtarget1
    real(rkind) :: dt1, dt2, dt, InflowSpeed = 1.d0
    real(rkind) :: x_shift
    integer :: nxADSIM, nxHIT
    type(fof), dimension(:), allocatable :: filt
    integer, dimension(:), allocatable :: pid
    integer :: fid, nfilters = 2, tid_FIL_FullField = 75, tid_FIL_Planes = 4
    logical :: applyFilters = .false. 

    namelist /concurrent/ XXXXXXXXXXXXXXXXXXXXXXXX

    call MPI_Init(ierr)                                                

    call GETARG(1,inputfile)                                            

    allocate(hit)                                                       
    allocate(adsim)                                                     
    
    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=concurrent)
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
    XXX need to deal with the temperature XXX             
    utarget0 = InflowSpeed                                                     
    vtarget0 = 0.d0                                                      
    wtarget0 = 0.d0                                                      
    call adsim%fringe_x1%associateFringeTargets(utarget0, vtarget0, wtarget0) 

    call adsim%fringe_x2%allocateTargetArray_Cells(utarget1)                
    call adsim%fringe_x2%allocateTargetArray_Cells(vtarget1)                
    call adsim%fringe_x2%allocateTargetArray_Edges(wtarget1)  
    XXX need to deal with the temperature XXX              
    utarget1 = InflowSpeed
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
    uTarget1(nxADSIM-nxHIT+1:nxADSIM,:,:) = uTarget1(nxADSIM-nxHIT+1:nxADSIM,:,:) + InflowSpeed 
    vTarget1(nxADSIM-nxHIT+1:nxADSIM,:,:) = vTarget1(nxADSIM-nxHIT+1:nxADSIM,:,:)
    wTarget1(nxADSIM-nxHIT+1:nxADSIM,:,:) = wTarget1(nxADSIM-nxHIT+1:nxADSIM,:,:)

    call budg_tavg%init(AD_Inputfile, adsim)   !<-- Budget class initialization 

    call message("==========================================================")
    call message(0, "All memory allocated! Now running the simulation.")
    call tic() 
    do while (adsim%tsim < adsim%tstop) 
       dt1 = adsim%get_dt(recompute=.true.)
       dt2 = hit%get_dt(recompute=.true.)
       dt = min(dt1, dt2)
       
       call adsim%timeAdvance(dt)

       call hit%timeAdvance(dt)
       
       call budg_tavg%doBudgets()       !<--- perform budget related operations -----Question::Where should this be placed ??

       x_shift = adsim%tsim*InflowSpeed 
       call hit%spectC%bandpassFilter_and_phaseshift(hit%whatC , uTarget1(nxADSIM-nxHIT+1:nxADSIM,:,:), x_shift)
       call hit%interpolate_cellField_to_edgeField(uTarget1(nxADSIM-nxHIT+1:nxADSIM,:,:), wTarget1(nxADSIM-nxHIT+1:nxADSIM,:,:),0,0)
       call hit%spectC%bandpassFilter_and_phaseshift(hit%uhat  , uTarget1(nxADSIM-nxHIT+1:nxADSIM,:,:), x_shift)
       call hit%spectC%bandpassFilter_and_phaseshift(hit%vhat  , vTarget1(nxADSIM-nxHIT+1:nxADSIM,:,:), x_shift)
      
       ! Now scale rhw HIT field appropriately
       uTarget1(nxADSIM-nxHIT+1:nxADSIM,:,:) = uTarget1(nxADSIM-nxHIT+1:nxADSIM,:,:) + InflowSpeed 
       vTarget1(nxADSIM-nxHIT+1:nxADSIM,:,:) = vTarget1(nxADSIM-nxHIT+1:nxADSIM,:,:)
       wTarget1(nxADSIM-nxHIT+1:nxADSIM,:,:) = wTarget1(nxADSIM-nxHIT+1:nxADSIM,:,:)
     

       call doTemporalStuff(adsim, 1)                                        
       call doTemporalStuff(hit  , 2)                                        
     
    end do 

    call budg_tavg%destroy()           !<-- release memory taken by the budget class 

    call hit%finalize_io()          
    call adsim%finalize_io()         

    call hit%destroy()                
    call adsim%destroy()               
   
    deallocate(hit, adsim)             
    
    call MPI_Finalize(ierr)        

end program
