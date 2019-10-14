subroutine instrumentForBudgets(this, uc, vc, wc, usgs, vsgs, wsgs, uvisc, vvisc, wvisc, px, py, pz, wb, ucor, vcor, wcor, uturb)  
    class(igrid), intent(inout) :: this
    complex(rkind), dimension(:,:,:), intent(in), target :: uc, vc, wc, usgs, vsgs, wsgs, uvisc, vvisc, wvisc, px, py, pz, wb, ucor, vcor, wcor, uturb 

    this%ucon => uc
    this%vcon => vc
    this%wcon => wc

    this%usgs => usgs
    this%vsgs => vsgs
    this%wsgs => wsgs

    this%uvisc => uvisc
    this%vvisc => vvisc
    this%wvisc => wvisc


    this%ucor => ucor
    this%vcor => vcor
    this%wcor => wcor

    this%px => px
    this%py => py
    this%pz => pz

    this%wb => wb
    this%uturb => uturb 

    ! Safeguards
    this%StoreForBudgets = .true. 
    if (.not. this%fastCalcPressure) then
        call GracefulExit("Cannot perform budget calculations if IGRID is initialized with FASTCALCPRESSURE=.false.", 324)
    end if 
  
    if (.not. useSkewSymm) then
        call message("WARNING: Advection term should be evaluated in the skew-symmetric form in order to perform budget calculations.")
    end if 

    if (this%useControl) then
        call message("WARNING: Budget calculations ignore the frame angle controller effects.", 324)
    end if 
    call message(0, "Budget calculations instrumented within igrid!")

    call this%set_budget_rhs_to_zero()
end subroutine

subroutine instrumentForBudgets_TimeAvg(this, uc, vc, wc, usgs, vsgs, wsgs,  px, py, pz, uturb, vturb, wturb, pxdns, pydns, pzdns, uvisc, vvisc, wvisc, ucor, vcor, wcor, wb)  
    class(igrid), intent(inout) :: this
    complex(rkind), dimension(:,:,:), intent(in), target :: uc, vc, wc, usgs, vsgs, wsgs, px, py, pz, uturb, vturb, wturb
    complex(rkind), dimension(:,:,:), intent(in), target :: pxdns, pydns, pzdns, uvisc, vvisc, wvisc, ucor, vcor, wcor, wb

    this%ucon => uc
    this%vcon => vc
    this%wcon => wc

    this%usgs => usgs
    this%vsgs => vsgs
    this%wsgs => wsgs

    this%px => px
    this%py => py
    this%pz => pz

    this%uvisc => uvisc
    this%vvisc => vvisc
    this%wvisc => wvisc
    
    this%ucor => ucor
    this%vcor => vcor
    this%wcor => wcor
    
    this%wb => wb

    this%pxdns => pxdns
    this%pydns => pydns
    this%pzdns => pzdns

    this%uturb => uturb 
    this%vturb => vturb 
    this%wturb => wturb 


    ! Safeguards
    this%StoreForBudgets = .true. 
    if (.not. this%fastCalcPressure) then
        call GracefulExit("Cannot perform budget calculations if IGRID is initialized with FASTCALCPRESSURE=.false.", 324)
    end if 
  
    if (.not. useSkewSymm) then
        call message("WARNING: Advection term should be evaluated in the skew-symmetric form in order to perform budget calculations.")
    end if 

    if (this%useControl) then
        call message("WARNING: Budget calculations ignore the frame angle controller effects.", 324)
    end if 
    call message(0, "Time averaging budget calculations instrumented within igrid!")

    call this%set_budget_rhs_to_zero()
end subroutine

subroutine set_budget_rhs_to_zero(this)
    use constants, only: imi

    class(igrid), intent(inout) :: this
    complex(rkind) :: czero 

    czero = 0.d0 + imi*0.d0
    
    this%ucon = czero 
    this%vcon = czero 
    this%wcon = czero 

    this%usgs = czero 
    this%vsgs = czero 
    this%wsgs = czero 

    this%px = czero 
    this%py = czero 
    this%pz = czero

    this%uvisc = czero 
    this%vvisc = czero 
    this%wvisc = czero 

    if (associated(this%ucor)) this%ucor = czero 
    if (associated(this%vcor)) this%vcor = czero 
    if (associated(this%wcor)) this%wcor = czero 

    if (associated(this%wb   )) this%wb = czero 
    if (associated(this%uturb)) this%uturb = czero 
    if (associated(this%vturb)) this%vturb = czero 
    if (associated(this%wturb)) this%wturb = czero 

end subroutine 

subroutine getMomentumTerms(this)
    class(igrid), intent(inout) :: this

    if (this%StoreForBudgets) then
        ! Step 1: split the RHS terms into 4 contributions and pressure
        
        
        call this%ComputePressure(ComputeRHSForBudget=.true.)
        
        ! Step 2: compute the pressure gradient term 
        call this%spectC%fft(this%Pressure, this%cbuffyC(:,:,:,1))
        this%cbuffyC(:,:,:,1) = -this%cbuffyC(:,:,:,1) ! Pressure terms is -gradP
        call this%spectC%mtimes_ik1_oop(this%cbuffyC(:,:,:,1),this%px)
        call this%spectC%mtimes_ik2_oop(this%cbuffyC(:,:,:,1),this%py)
        call transpose_y_to_z(this%cbuffyC(:,:,:,1),this%cbuffzC(:,:,:,1),this%sp_gpC)
        call this%Pade6opZ%ddz_C2E(this%cbuffzC(:,:,:,1),this%cbuffzE(:,:,:,1),0,0)  ! Safest to do 0,0 for BC because it's most general
        this%cbuffyE(:,:,:,1) = this%wcon + this%wsgs
        if(associated(this%wb)) this%cbuffyE(:,:,:,1) = this%cbuffyE(:,:,:,1) + this%wb
        if(associated(this%wcor)) this%cbuffyE(:,:,:,1) = this%cbuffyE(:,:,:,1) + this%wcor
        ! Set the true BC for dpdz:
        call transpose_y_to_z(this%cbuffyE(:,:,:,1),this%cbuffzE(:,:,:,2),this%sp_gpE)
        this%cbuffzE(:,:,1,1) = -this%cbuffzE(:,:,1,2)
        this%cbuffzE(:,:,this%nz+1,1) = -this%cbuffzE(:,:,this%nz+1,2)
        call transpose_z_to_y(this%cbuffzE(:,:,:,1),this%pz, this%sp_gpE)

        if(associated(this%pxdns) .and. associated(this%pydns) .and. associated(this%pzdns)) then
            call this%spectC%fft(this%pressure_dns, this%cbuffyC(:,:,:,1))
            this%cbuffyC(:,:,:,1) = -this%cbuffyC(:,:,:,1) ! Pressure terms is -gradP
            call this%spectC%mtimes_ik1_oop(this%cbuffyC(:,:,:,1),this%pxdns)
            call this%spectC%mtimes_ik2_oop(this%cbuffyC(:,:,:,1),this%pydns)
            call transpose_y_to_z(this%cbuffyC(:,:,:,1),this%cbuffzC(:,:,:,1),this%sp_gpC)
            call this%Pade6opZ%ddz_C2E(this%cbuffzC(:,:,:,1),this%cbuffzE(:,:,:,1),0,0)  ! Safest to do 0,0 for BC because it's most general
            call transpose_z_to_y(this%cbuffzE(:,:,:,1),this%pzdns, this%sp_gpE)
        endif

    end if 

end subroutine 
