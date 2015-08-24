#include "hitSpectral_files/variables.F90"
#include "hitSpectral_files/io.F90"
#include "hitSpectral_files/spectralOps.F90"
#include "hitSpectral_files/initialization.F90"
#include "hitSpectral_files/timestep.F90"

program hitSpectral
    use kind_parameters, only: rkind
    use initialization,  only: initialize, finalize  
    use variables,       only: maxDivergence, NT, FT, fieldsSpec, fieldsPhys, RESTART, dt, rhs,time,&
                                    BUFF_REAL, TSTEP_DUMP, TSTEP_RESTART, rhsOld, DealiasMat

    use decomp_2d,       only: nrank  
    use timeStep,        only: getDT, WrapUpTstep 
    use io,              only: dumpData4Matlab, createRestartFile
    use spectralOps,     only: divergence, PressureProjection, getrhs 

    implicit none
    integer :: tidx, tidx_start
    logical :: checkDivergence = .false.     

    ! Initialize all fields (read input files)
    call initialize(tidx_start) 

    ! Get time step from CFL (or use DT) 
    ! NOTE: spectral -> physical transform occurs in this step
    call getDT

    
    ! First time step after initialization/restart 
    call getrhs
    if (.not. RESTART) then
        fieldsSpec = fieldsSpec + dt*rhs
        call pressureProjection
    else
        fieldsSpec = fieldsSpec + dt*(1.5_rkind*rhs - 0.5_rkind*rhsOld)
        call pressureProjection
    end if
    ! Dealias 
    fieldsSpec(:,:,:,1) = DealiasMat*fieldsSpec(:,:,:,1)
    fieldsSpec(:,:,:,2) = DealiasMat*fieldsSpec(:,:,:,2)
    fieldsSpec(:,:,:,3) = DealiasMat*fieldsSpec(:,:,:,3)
    

    ! Prepare for next time step 
    rhsOld = rhs
    time = time + dt
    call WrapUpTstep(tidx_start)! Convert u from Spec-> Phys 
    tidx_start = tidx_start + 1
    
    
    ! Update time step
    call getDT
    ! Regular time stepping using Adams Bashforth
    do tidx = tidx_start,NT
        ! Do the actual time advancement
        call getrhs
        fieldsSpec = fieldsSpec + dt*(1.5_rkind*rhs - 0.5_rkind*rhsOld)
        call pressureProjection
       
        ! OPTIONAL - Check if the field is solenoidal  
        if (checkDivergence) then
            call divergence
            if(nrank == 0) print*, maxDivergence
        end if

        ! DEALIAS in preparation for the next step  
        fieldsSpec(:,:,:,1) = DealiasMat*fieldsSpec(:,:,:,1)
        fieldsSpec(:,:,:,2) = DealiasMat*fieldsSpec(:,:,:,2)
        fieldsSpec(:,:,:,3) = DealiasMat*fieldsSpec(:,:,:,3)
    
        ! Prepare for next time step  
        RHSold = RHS
        time = time + dt 
        call getDT 
        call WrapUpTstep(tidx) 
        
        
        ! File IO operations
        if (mod(tidx,TSTEP_DUMP) == 0) then
            call dumpData4Matlab(tidx)
        end if 
        if (mod(tidx,TSTEP_RESTART) == 0) then
           call createRestartFile 
        end if 

    end do 

    call finalize 
end program hitSpectral
