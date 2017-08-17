! Template for PadeOps

#include "HIT_AD_interact_files/initialize.F90"       
#include "HIT_AD_interact_files/temporalHook.F90"  

program HIT_AD_interact
    use mpi
    use kind_parameters,  only: clen, rkind
    use IncompressibleGrid, only: igrid
    use temporalhook, only: doTemporalStuff
    use timer, only: tic, toc
    use exits, only: message
    use constants, only: one, zero
    use HIT_AD_interact_parameters, only: simulationID
    use decomp_2d,                only: nrank
    implicit none

    type(igrid), allocatable, target :: hit, adsim
    character(len=clen) :: inputfile, HIT_InputFile, AD_InputFile
    integer :: ierr, ioUnit
    real(rkind), dimension(:,:,:), allocatable :: utarget0, vtarget0, wtarget0
    real(rkind), dimension(:,:,:), allocatable :: utarget1, vtarget1, wtarget1
    real(rkind) :: HIT_scalefact = 0.2d0, HIT_TI = 0.1d0, dt1, dt2, dt
    real(rkind) :: k_bandpass_left = 10.d0, k_bandpass_right = 64.d0, x_shift
    integer :: nxADSIM, nxHIT
    namelist /concurrent/ HIT_InputFile, AD_InputFile, HIT_scalefact, HIT_TI, k_bandpass_left, k_bandpass_right

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
    call hit%start_io(.true.)                                           
    call hit%printDivergence()
    
    call adsim%fringe_x1%allocateTargetArray_Cells(utarget0)                
    call adsim%fringe_x1%allocateTargetArray_Cells(vtarget0)                
    call adsim%fringe_x1%allocateTargetArray_Edges(wtarget0)                
    utarget0 = 1.d0                                                      
    vtarget0 = 0.d0                                                      
    wtarget0 = 0.d0                                                      
    call adsim%fringe_x1%associateFringeTargets(utarget0, vtarget0, wtarget0) 

    call adsim%fringe_x2%allocateTargetArray_Cells(utarget1)                
    call adsim%fringe_x2%allocateTargetArray_Cells(vtarget1)                
    call adsim%fringe_x2%allocateTargetArray_Edges(wtarget1)                
    utarget1 = 1.d0
    vtarget1 = 0.d0
    wtarget1 = 0.d0
    call adsim%fringe_x2%associateFringeTargets(utarget1, vtarget1, wtarget1) 

    nxHIT = hit%gpC%xsz(1)
    nxADSIM  = adsim%gpC%xsz(1)

    call hit%spectC%init_bandpass_filter(k_bandpass_left, k_bandpass_right, hit%cbuffzC(:,:,:,1), hit%cbuffyC(:,:,:,1))

    ! Set the true target field for AD simulation 
    call hit%spectC%bandpassFilter_and_phaseshift(hit%whatC , uTarget1(nxADSIM-nxHIT+1:nxADSIM,:,:), zero)
    call hit%interpolate_cellField_to_edgeField(uTarget1(nxADSIM-nxHIT+1:nxADSIM,:,:), wTarget1(nxADSIM-nxHIT+1:nxADSIM,:,:),0,0)
    call hit%spectC%bandpassFilter_and_phaseshift(hit%uhat  , uTarget1(nxADSIM-nxHIT+1:nxADSIM,:,:), zero)
    call hit%spectC%bandpassFilter_and_phaseshift(hit%vhat  , vTarget1(nxADSIM-nxHIT+1:nxADSIM,:,:), zero)
    
    ! Now scale rhw HIT field appropriately
    ! Note that the bandpass filtered velocity field has zero mean 
    uTarget1(nxADSIM-nxHIT+1:nxADSIM,:,:) = HIT_scalefact*uTarget1(nxADSIM-nxHIT+1:nxADSIM,:,:) + one 
    vTarget1(nxADSIM-nxHIT+1:nxADSIM,:,:) = HIT_scalefact*vTarget1(nxADSIM-nxHIT+1:nxADSIM,:,:)
    wTarget1(nxADSIM-nxHIT+1:nxADSIM,:,:) = HIT_scalefact*wTarget1(nxADSIM-nxHIT+1:nxADSIM,:,:)
    
    call message("==========================================================")
    call message(0, "All memory allocated! Now running the simulation.")
    call tic() 
    do while (adsim%tsim < adsim%tstop) 
       dt1 = adsim%get_dt(recompute=.true.)
       dt2 = hit%get_dt(recompute=.true.)
       dt = min(dt1, dt2)
       
       call adsim%timeAdvance(dt)

       call hit%timeAdvance(dt)
       
       call hit%spectC%bandpassFilter_and_phaseshift(hit%whatC , uTarget1(nxADSIM-nxHIT+1:nxADSIM,:,:), zero)
       call hit%interpolate_cellField_to_edgeField(uTarget1(nxADSIM-nxHIT+1:nxADSIM,:,:), wTarget1(nxADSIM-nxHIT+1:nxADSIM,:,:),0,0)
       call hit%spectC%bandpassFilter_and_phaseshift(hit%uhat  , uTarget1(nxADSIM-nxHIT+1:nxADSIM,:,:), zero)
       call hit%spectC%bandpassFilter_and_phaseshift(hit%vhat  , vTarget1(nxADSIM-nxHIT+1:nxADSIM,:,:), zero)
      
       ! Now scale rhw HIT field appropriately
       uTarget1(nxADSIM-nxHIT+1:nxADSIM,:,:) = HIT_scalefact*uTarget1(nxADSIM-nxHIT+1:nxADSIM,:,:) + one 
       vTarget1(nxADSIM-nxHIT+1:nxADSIM,:,:) = HIT_scalefact*vTarget1(nxADSIM-nxHIT+1:nxADSIM,:,:)
       wTarget1(nxADSIM-nxHIT+1:nxADSIM,:,:) = HIT_scalefact*wTarget1(nxADSIM-nxHIT+1:nxADSIM,:,:)
     

       call doTemporalStuff(adsim, 1)                                        
       call doTemporalStuff(hit   , 2)                                        
      
    end do 
 
    call hit%finalize_io()          
    call adsim%finalize_io()         

    call hit%destroy()                
    call adsim%destroy()               
   
    deallocate(hit, adsim)             
    
    call MPI_Finalize(ierr)        

end program
