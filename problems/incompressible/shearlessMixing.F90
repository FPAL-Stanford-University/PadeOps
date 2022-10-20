! Template for PadeOps

#include "shearlessMixing_files/initializeV2.F90"       
#include "shearlessMixing_files/temporalHook.F90"  

program shearlessMixing
    use mpi
    use kind_parameters,                     only: clen, rkind
    use IncompressibleGrid,                  only: igrid
    use temporalhook,                        only: doTemporalStuff
    use timer,                               only: tic, toc
    use exits,                               only: message, GracefulExit
    use fortran_assert,                      only: assert
    use constants,                           only: one, zero
    use shearlessMixing_interact_parameters, only: simulationID, &
                                             streamWiseCoord, nxHIT, nyHIT, nzHIT, &
                                             nxSM, nySM, nzSM, copyHITfieldsToSM, &
                                             MMS, finalizeSimulation, checkMMS
    use fof_mod,                             only: fof

    implicit none

    type(igrid), allocatable, target :: HIT, SM
    character(len=clen) :: inputfile, HIT_InputFile, SM_InputFile, fof_dir, filoutdir
    integer :: ierr, ioUnit

    real(rkind), dimension(:,:,:), allocatable :: utarget0, vtarget0, wtarget0, Ttarget0
    real(rkind), dimension(:,:,:), allocatable :: uHITFilt, vHITFilt, wHITFilt
    
    real(rkind) :: dt1, dt2, dt
    real(rkind) :: k_bandpass_left = 10.d0, k_bandpass_right = 64.d0
    type(fof), dimension(:), allocatable :: filt
    integer, dimension(:), allocatable :: pid
    integer :: fid, nfilters = 2, tid_FIL_FullField = 75, tid_FIL_Planes = 4
    logical :: applyFilters = .false. 
    logical, parameter :: synchronize_RK_substeps = .true.

    namelist /concurrent/ HIT_InputFile, SM_InputFile, k_bandpass_left, &
      k_bandpass_right, MMS

    call MPI_Init(ierr)                                                

    call GETARG(1,inputfile)                                            

    allocate(HIT)                                                       
    allocate(SM)                                                     
    
    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=concurrent)
    close(ioUnit)    

    simulationID = 1
    call SM%init(SM_InputFile, .true.)                               
    call SM%start_io(.true.)                                          
    call SM%printDivergence()
   
    call mpi_barrier(mpi_comm_world, ierr)

    simulationID = 2
    call HIT%init(HIT_InputFile, .false.)                                           
    HIT%Am_I_Primary = .false. 
    call HIT%start_io(.true.)                                           
    call HIT%printDivergence()

    allocate(uHITFilt(HIT%gpC%xsz(1),HIT%gpC%xsz(2),HIT%gpC%xsz(3)))
    allocate(vHITFilt(HIT%gpC%xsz(1),HIT%gpC%xsz(2),HIT%gpC%xsz(3)))
    allocate(wHITFilt(HIT%gpE%xsz(1),HIT%gpE%xsz(2),HIT%gpE%xsz(3)))
    
    call SM%fringe_x%allocateTargetArray_Cells(utarget0)                
    call SM%fringe_x%allocateTargetArray_Cells(vtarget0)                
    call SM%fringe_x%allocateTargetArray_Cells(Ttarget0)                
    call SM%fringe_x%allocateTargetArray_Edges(wtarget0)                
    call SM%fringe_x%associateFringeTargets(utarget0, vtarget0, wtarget0) 
    call SM%fringe_x%associateFringeTarget_scalar(Ttarget0) 


    nxHIT = HIT%gpC%xsz(1)
    nyHIT = HIT%gpC%ysz(2)
    nzHIT = HIT%gpC%zsz(3)
    nxSM  = SM%gpC%xsz(1)
    nySM  = SM%gpC%ysz(2)
    nzSM  = SM%gpC%zsz(3)

    call HIT%spectC%init_bandpass_filter(k_bandpass_left, k_bandpass_right, HIT%cbuffzC(:,:,:,1), HIT%cbuffyC(:,:,:,1))

    ! Set the true target field for AD simulation
    call HIT%spectC%bandpassFilter_and_phaseshift(HIT%whatC , uHITFilt, zero, streamWiseCoord)
    call HIT%interpolate_cellField_to_edgeField(uHITFilt,wHITFilt,0,0)
    call HIT%spectC%bandpassFilter_and_phaseshift(HIT%uhat  , uHITFilt, zero, streamWiseCoord)
    call HIT%spectC%bandpassFilter_and_phaseshift(HIT%vhat  , vHITFilt, zero, streamWiseCoord)
    call copyHITfieldsToSM(uHITFilt,vHITFilt,wHITFilt,utarget0,vtarget0,wtarget0,HIT,SM,streamWiseCoord)

    ! Initialize the filters
    if (applyfilters) then
      allocate(filt(nfilters))
      do fid = 1,nfilters
         call filt(fid)%init(SM%runID, fof_dir, filoutdir, fid, SM%spectC, SM%cbuffyC(:,:,:,1), SM%cbuffzC(:,:,:,1), .true., SM%gpC)
      end do 
      if (allocated(SM%yplanes)) then 
         allocate(pid(size(SM%yplanes)))
         pid = SM%yplanes 
      end if 
    end if 


    call message("==========================================================")
    call message(0, "All memory allocated! Now running the simulation.")
    call tic() 
    do while (SM%tsim < SM%tstop) 
       dt1 = SM%get_dt(recompute=.true.)
       dt2 = HIT%get_dt(recompute=.true.)
       dt = min(dt1, dt2)
      
       SM%dt = dt
       HIT%dt = dt

       select case (HIT%TimeSteppingScheme)
       case (0)
         call assert(SM%TimeSteppingScheme == 0, 'SM%TimeSteppingScheme == 0')
         call SM%timeAdvance(dt)
         call HIT%timeAdvance(dt)

         call HIT%spectC%bandpassFilter_and_phaseshift(HIT%whatC, uHITFilt, &
           zero, streamWiseCoord)
         call HIT%interpolate_cellField_to_edgeField(uHITFilt, wHITFilt,0,0)
         call HIT%spectC%bandpassFilter_and_phaseshift(HIT%uhat, uHITFilt, &
           zero, streamWiseCoord)
         call HIT%spectC%bandpassFilter_and_phaseshift(HIT%vhat, vHITFilt, &
           zero, streamWiseCoord)
         call copyHITfieldsToSM(uHITFilt,vHITFilt,wHITFilt,utarget0,vtarget0,&
           wtarget0,HIT,SM,streamWiseCoord)
       case (2)
         ! Stage 1
         call SM%advance_SSP_RK45_Stage_1()
         call HIT%advance_SSP_RK45_Stage_1()
         call HIT%spectC%bandpassFilter_and_phaseshift(HIT%whatC, uHITFilt, &
           zero, streamWiseCoord)
         call HIT%interpolate_cellField_to_edgeField(uHITFilt, wHITFilt,0,0)
         call HIT%spectC%bandpassFilter_and_phaseshift(HIT%uhat, uHITFilt, &
           zero, streamWiseCoord)
         call HIT%spectC%bandpassFilter_and_phaseshift(HIT%vhat, vHITFilt, &
           zero, streamWiseCoord)
         call copyHITfieldsToSM(uHITFilt,vHITFilt,wHITFilt,utarget0,vtarget0,&
           wtarget0,HIT,SM,streamWiseCoord)

         ! Stage 2
         call SM%advance_SSP_RK45_Stage_2()
         call HIT%advance_SSP_RK45_Stage_2()
         call HIT%spectC%bandpassFilter_and_phaseshift(HIT%whatC, uHITFilt, &
           zero, streamWiseCoord)
         call HIT%interpolate_cellField_to_edgeField(uHITFilt, wHITFilt,0,0)
         call HIT%spectC%bandpassFilter_and_phaseshift(HIT%uhat, uHITFilt, &
           zero, streamWiseCoord)
         call HIT%spectC%bandpassFilter_and_phaseshift(HIT%vhat, vHITFilt, &
           zero, streamWiseCoord)
         call copyHITfieldsToSM(uHITFilt,vHITFilt,wHITFilt,utarget0,vtarget0,&
           wtarget0,HIT,SM,streamWiseCoord)

         ! Stage 3
         call SM%advance_SSP_RK45_Stage_3()
         call HIT%advance_SSP_RK45_Stage_3()
         call HIT%spectC%bandpassFilter_and_phaseshift(HIT%whatC, uHITFilt, &
           zero, streamWiseCoord)
         call HIT%interpolate_cellField_to_edgeField(uHITFilt, wHITFilt,0,0)
         call HIT%spectC%bandpassFilter_and_phaseshift(HIT%uhat, uHITFilt, &
           zero, streamWiseCoord)
         call HIT%spectC%bandpassFilter_and_phaseshift(HIT%vhat, vHITFilt, &
           zero, streamWiseCoord)
         call copyHITfieldsToSM(uHITFilt,vHITFilt,wHITFilt,utarget0,vtarget0,&
           wtarget0,HIT,SM,streamWiseCoord)

         ! Stage 4
         call SM%advance_SSP_RK45_Stage_4()
         call HIT%advance_SSP_RK45_Stage_4()
         call HIT%spectC%bandpassFilter_and_phaseshift(HIT%whatC, uHITFilt, &
           zero, streamWiseCoord)
         call HIT%interpolate_cellField_to_edgeField(uHITFilt, wHITFilt,0,0)
         call HIT%spectC%bandpassFilter_and_phaseshift(HIT%uhat, uHITFilt, &
           zero, streamWiseCoord)
         call HIT%spectC%bandpassFilter_and_phaseshift(HIT%vhat, vHITFilt, &
           zero, streamWiseCoord)
         call copyHITfieldsToSM(uHITFilt,vHITFilt,wHITFilt,utarget0,vtarget0,&
           wtarget0,HIT,SM,streamWiseCoord)

         ! Stage 5
         call SM%advance_SSP_RK45_Stage_5()
         call HIT%advance_SSP_RK45_Stage_5()
         call HIT%spectC%bandpassFilter_and_phaseshift(HIT%whatC, uHITFilt, &
           zero, streamWiseCoord)
         call HIT%interpolate_cellField_to_edgeField(uHITFilt, wHITFilt,0,0)
         call HIT%spectC%bandpassFilter_and_phaseshift(HIT%uhat, uHITFilt, &
           zero, streamWiseCoord)
         call HIT%spectC%bandpassFilter_and_phaseshift(HIT%vhat, vHITFilt, &
           zero, streamWiseCoord)
         call copyHITfieldsToSM(uHITFilt,vHITFilt,wHITFilt,utarget0,vtarget0,&
           wtarget0,HIT,SM,streamWiseCoord)

       ! Call wrap up 
         call SM%wrapup_timestep()
         call HIT%wrapup_timestep() 
       case default
         call gracefulExit('Time stepping scheme not supported for this problem',ierr)
       end select

       if (MMS) then
         call checkMMS(SM%u, SM%v, SM%w, SM%mesh, SM%tsim)
       end if

       call doTemporalStuff(SM, 1)                                        
       call doTemporalStuff(HIT  , 2)                                        
     
    end do

    deallocate(uHITFilt, vHITFilt, wHITFilt)
    deallocate(utarget0, vtarget0, wtarget0)
    call HIT%finalize_io()
    call SM%finalize_io()

    call HIT%destroy()
    call SM%destroy()
   
    deallocate(HIT, SM)

    call finalizeSimulation
    
    call MPI_Finalize(ierr)

end program
