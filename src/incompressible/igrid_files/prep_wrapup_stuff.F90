subroutine dealiasFields(this)
    class(igrid), intent(inout) :: this
    integer :: idx

    if (this%donot_dealias) then
       return
    else
        call this%spectC%dealias(this%uhat)
        call this%spectC%dealias(this%vhat)
        if (this%PeriodicInZ) then
            call transpose_y_to_z(this%what, this%cbuffzE(:,:,:,1), this%sp_gpE)
            call this%spectC%dealias_edgeField(this%cbuffzE(:,:,:,1))
            call transpose_z_to_y(this%cbuffzE(:,:,:,1),this%what,this%sp_gpE)
        else
            call this%spectE%dealias(this%what)
        end if 
        if (this%isStratified .or. this%initspinup) call this%spectC%dealias(this%That)
        if (this%UseDealiasFilterVert) then
            call this%ApplyCompactFilter()
        end if

        if ((this%usescalars) .and. allocated(this%scalars)) then
           do idx = 1,this%n_scalars
              call this%scalars(idx)%dealias()
           end do 
        end if 
    end if 
end subroutine 


 subroutine dealias_rhs(this, uin, vin, win)
    class(igrid), intent(inout) :: this
    complex(rkind), dimension(this%sp_gpC%ysz(1),this%sp_gpC%ysz(2),this%sp_gpC%ysz(3)), intent(inout ) :: uin, vin
    complex(rkind), dimension(this%sp_gpE%ysz(1),this%sp_gpE%ysz(2),this%sp_gpE%ysz(3)), intent(inout) :: win

   if (this%donot_dealias) then
      return
   else
      call this%spectC%dealias(uin)
      call this%spectC%dealias(vin)
      if (this%PeriodicInZ) then
          call transpose_y_to_z(win, this%cbuffzE(:,:,:,1), this%sp_gpE)
          call this%spectC%dealias_edgeField(this%cbuffzE(:,:,:,1))
          call transpose_z_to_y(this%cbuffzE(:,:,:,1),win,this%sp_gpE)
      else
          call this%spectE%dealias(win)
      end if 
   end if 
end subroutine 

subroutine dealiasRealField_C(this, field)
    class(igrid), intent(inout) :: this
    real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(1),this%gpC%xsz(1)), intent(inout) :: field
  
    if (this%donot_dealias) then
        return 
    else
        call this%spectC%fft(field, this%cbuffyC(:,:,:,1))
        call this%spectC%dealias(this%cbuffyC(:,:,:,1))
        call this%spectC%ifft(this%cbuffyC(:,:,:,1),field)
    end if 
end subroutine

subroutine computePressure(this, ComputeRHSForBudget)
    class(igrid), intent(inout) :: this
    logical :: copyFringeRHS, copyDNSRHS, copyTurbineRHS
    logical, optional, intent(in) :: ComputeRHSForBudget
    logical :: splitRHS

    if (present(ComputeRHSForBudget)) then
        splitRHS = ComputeRHSForBudget
    else
        splitRHS = .false. 
    end if 

    if (this%computeDNSpressure) then
       copyDNSRHS = .true. 
    else
       copyDNSRHS = .false. 
    end if

    if (this%computeFringePressure) then
       copyFringeRHS = .true. 
    else
       copyFringeRHS = .false. 
    end if

    if (this%computeTurbinePressure) then
        copyTurbineRHS = .true.
     else
        copyTurbineRHS = .false. 
    end if

    ! STEP 1: Populate RHS
    if (splitRHS) then
        call this%populate_rhs_for_budgets(CopyForDNSpress=copyDNSRHS, CopyForFringePress=copyFringeRHS, copyForTurbinePress=copyTurbineRHS)
    else
        call this%populate_rhs(CopyForDNSpress=copyDNSRHS, CopyForFringePress=copyFringeRHS, copyForTurbinePress=copyTurbineRHS)
    end if

    ! STEP 2: Compute pressure
    if (this%fastCalcPressure) then
       call this%dealias_rhs(this%u_rhs, this%v_rhs, this%w_rhs) 
       call this%Padepoiss%getPressureAndUpdateRHS(this%u_rhs,this%v_rhs,this%w_rhs,this%pressure)
        !call this%dealiasRealField_C(this%pressure)
        if (this%computeDNSpressure) then
           call this%dealias_rhs(this%urhs_dns, this%vrhs_dns, this%wrhs_dns) 
           call this%padepoiss%getPressure(this%urhs_dns,this%vrhs_dns,this%wrhs_dns,this%pressure_dns)
           !call this%dealiasRealField_C(this%pressure_dns)
        end if
        if (this%computeFringePressure) then
           call this%dealias_rhs(this%urhs_fringe, this%vrhs_fringe, this%wrhs_fringe) 
           call this%padepoiss%getPressure(this%urhs_fringe,this%vrhs_fringe,this%wrhs_fringe,this%pressure_fringe)
           !call this%dealiasRealField_C(this%pressure_fringe)
        end if
        if (this%computeTurbinePressure) then
           call this%dealias_rhs(this%urhs_turbine, this%vrhs_turbine, this%wrhs_turbine) 
           call this%padepoiss%getPressure(this%urhs_turbine,this%vrhs_turbine,this%wrhs_turbine,this%pressure_turbine)
           !call this%dealiasRealField_C(this%pressure_turbine)
        end if
    else
        call this%dealias_rhs(this%u_rhs, this%v_rhs, this%w_rhs) 
        call this%padepoiss%getPressure(this%u_rhs,this%v_rhs,this%w_rhs,this%pressure)
        !call this%dealiasRealField_C(this%pressure)
    end if 

    if (.not. useSkewSymm) then 
        ! You are using the rotational form. 
        ! This means that the pressure is really 
        ! the Bernoulli pressure. Need to subtract 
        ! out the kinetic energy. 
        call this%correctPressureRotationalForm()
    end if

    if (this%computeRapidSlowPressure) then
        call this%compute_RapidSlowPressure_Split()
    end if 
    
    ! STEP 3: Inform the other subroutines that you already have RHS
    this%AlreadyHaveRHS = .true. 

end subroutine

function getMaxKE(this) result(maxKE)
    class(igrid), intent(inout) :: this
    real(rkind)  :: maxKE

    this%rbuffxC(:,:,:,1) = this%u**2 + this%v**2 + this%wC**2
    maxKE = half*p_maxval(maxval(this%rbuffxC(:,:,:,1)))

end function

function getMeanKE(this) result(meanTKE)
     class(igrid), intent(inout) :: this
     real(rkind)  :: meanTKE

     this%rbuffxC(:,:,:,1) = this%u**2 + this%v**2 + this%wC**2
     meanTKE = half*p_sum(sum(this%rbuffxC(:,:,:,1)))/(real(this%nx)*real(this%ny)*real(this%nz))

end function

subroutine updateProbes(this)
    class(igrid), intent(inout) :: this
    integer :: idx

    if (this%doIhaveAnyProbes) then
        do idx = 1,this%nprobes
            this%probe_data(1,idx,this%step) = this%tsim
            this%probe_data(2,idx,this%step) = this%u (this%probes(1,idx),this%probes(2,idx),this%probes(3,idx))
            this%probe_data(3,idx,this%step) = this%v (this%probes(1,idx),this%probes(2,idx),this%probes(3,idx))
            this%probe_data(4,idx,this%step) = this%wC(this%probes(1,idx),this%probes(2,idx),this%probes(3,idx))
            if (this%isStratified) then
                this%probe_data(5,idx,this%step) = this%T(this%probes(1,idx),this%probes(2,idx),this%probes(3,idx))
            end if
            if (this%fastCalcPressure) then
                this%probe_data(6,idx,this%step) = this%Pressure(this%probes(1,idx),this%probes(2,idx),this%probes(3,idx))
            end if
            if (this%computeDNSpressure) then
                this%probe_data(7,idx,this%step) = this%Pressure_dns(this%probes(1,idx),this%probes(2,idx),this%probes(3,idx))
            end if 

            if (this%computeFringePressure) then
                this%probe_data(8,idx,this%step) = this%Pressure_fringe(this%probes(1,idx),this%probes(2,idx),this%probes(3,idx))
            end if
            
            if (this%computeTurbinePressure) then
                this%probe_data(9,idx,this%step) = this%Pressure_turbine(this%probes(1,idx),this%probes(2,idx),this%probes(3,idx))
            end if
        end do 
    end if

    ! KS - preprocess
    if (this%PreprocessForKS) then
        if (.not. this%KSupdated) then
            call this%LES2KS%applyFilterForKS(this%u, this%v, this%wC)
            this%KSupdated = .true. 
        end if
        if (this%doIhaveAnyProbes) then
            do idx = 1,this%nprobes
                this%KS_probe_data(1,idx,this%step) = this%tsim
                this%KS_probe_data(2,idx,this%step) = this%ufil4KS (this%probes(1,idx),this%probes(2,idx),this%probes(3,idx))
                this%KS_probe_data(3,idx,this%step) = this%vfil4KS (this%probes(1,idx),this%probes(2,idx),this%probes(3,idx))
                this%KS_probe_data(4,idx,this%step) = this%wfil4KS(this%probes(1,idx),this%probes(2,idx),this%probes(3,idx))
            end do 
        end if
    end if

end subroutine

subroutine correctPressureRotationalForm(this)
    class(igrid), intent(inout) :: this 
    real(rkind) :: mfact, meanK

    this%rbuffxC(:,:,:,1) = 0.5d0*(this%u*this%u + this%v*this%v + this%wC*this%wC)
    mfact = one/(real(this%nx,rkind)*real(this%ny,rkind)*real(this%nz,rkind))
    meanK = p_sum(sum(this%rbuffxC(:,:,:,1)))*mfact
    this%pressure = this%pressure + meanK
    this%pressure = this%pressure - this%rbuffxC(:,:,:,1)

    if (this%computeDNSpressure) then
        this%pressure_dns = this%pressure_dns + meanK
        this%pressure_dns = this%pressure_dns - this%rbuffxC(:,:,:,1)
    end if

end subroutine
   
subroutine interpolate_cellField_to_edgeField(this, rxC, rxE, bc1, bc2)
       class(igrid), intent(inout) :: this 
       real(rkind), intent(in),  dimension(this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3)) :: rxC 
       real(rkind), intent(out), dimension(this%gpE%xsz(1), this%gpE%xsz(2), this%gpE%xsz(3)) :: rxE
       integer, intent(in) :: bc1, bc2

       call transpose_x_to_y(rxC, this%rbuffyC(:,:,:,1), this%gpC)
       call transpose_y_to_z(this%rbuffyC(:,:,:,1), this%rbuffzC(:,:,:,1), this%gpC)
       call this%Pade6opZ%interpz_E2C(this%rbuffzC(:,:,:,1), this%rbuffzE(:,:,:,1), bc1,bc2) 
       call transpose_z_to_y(this%rbuffzE(:,:,:,1), this%rbuffyE(:,:,:,1), this%gpE)
       call transpose_y_to_x(this%rbuffyE(:,:,:,1), rxE, this%gpE)

end subroutine 

subroutine interp_PrimitiveVars(this)
    class(igrid), intent(inout), target :: this
    complex(rkind), dimension(:,:,:), pointer :: ybuffC, zbuffC, zbuffE
    
    zbuffE => this%cbuffzE(:,:,:,1)
    zbuffC => this%cbuffzC(:,:,:,1)
    ybuffC => this%cbuffyC(:,:,:,1)

    ! Step 1: Interpolate w -> wC
    call transpose_y_to_z(this%what,zbuffE,this%sp_gpE)
    call this%Pade6opZ%interpz_E2C(zbuffE,zbuffC,wBC_bottom, wBC_top)
    call transpose_z_to_y(zbuffC,this%whatC,this%sp_gpC)
    call this%spectC%ifft(this%whatC,this%wC)

    ! Step 2: Interpolate u -> uE
    call transpose_y_to_z(this%uhat,zbuffC,this%sp_gpC)
    call this%Pade6opZ%interpz_C2E(zbuffC,zbuffE,uBC_bottom, uBC_top)
    if (botWall==1) zbuffE(:,:,1)         = 0.d0 ! No-slip wall (sided scheme is corrected) 
    if (topWall==1) zbuffE(:,:,this%nz+1) = 0.d0 ! No-slip wall (sided scheme is corrected)  
    call transpose_z_to_y(zbuffE,this%uEhat, this%sp_gpE)
    call this%spectE%ifft(this%uEhat, this%uE)

    ! Step 3: Interpolate v -> vE
    call transpose_y_to_z(this%vhat,zbuffC,this%sp_gpC)
    call this%Pade6opZ%interpz_C2E(zbuffC,zbuffE,vBC_bottom, vBC_top)
    if (botWall==1) zbuffE(:,:,1)         = 0.d0 ! No-slip wall (sided scheme is corrected)
    if (topWall==1) zbuffE(:,:,this%nz+1) = 0.d0 ! No-slip wall (sided scheme is corrected) 
    call transpose_z_to_y(zbuffE,this%vEhat, this%sp_gpE)
    call this%spectE%ifft(this%vEhat, this%vE)
    

    ! Step 4: Interpolate T
    if (this%isStratified .or. this%initspinup) then
        call transpose_y_to_z(this%That,zbuffC,this%sp_gpC)
        call this%Pade6opZ%interpz_C2E(zbuffC,zbuffE,TBC_bottom, TBC_top)
        if ((this%botBC_Temp == 0)) then 
            zbuffE(:,:,1) = zero 
            if (nrank == 0) then
                zbuffE(1,1,1) = this%Tsurf*real(this%nx,rkind)*real(this%ny,rkind)
            end if 
        elseif ((this%botBC_Temp == 2)) then
            if (nrank == 0) then
                zbuffE(1,1,1) = this%Tsurf*real(this%nx,rkind)*real(this%ny,rkind)
            end if 
        end if

        ! Symmetric homogeneous dirichlet BC case
        if ((this%topBC_Temp == 3) .and. (this%botBC_Temp == 3)) then
            zbuffE(:,:,1) = zero
            zbuffE(:,:,this%nz+1) = zero 
        end if 
        call transpose_z_to_y(zbuffE,this%TEhat,this%sp_gpE)
        call this%spectE%ifft(this%TEhat,this%TE)
    end if 
end subroutine
   
subroutine ApplyCompactFilter(this)
    class(igrid), intent(inout), target :: this
    complex(rkind), dimension(:,:,:), pointer :: zbuff1, zbuff2, zbuff3, zbuff4
    zbuff1 => this%cbuffzC(:,:,:,1)
    zbuff2 => this%cbuffzC(:,:,:,2)
    zbuff3 => this%cbuffzE(:,:,:,1)
    zbuff4 => this%cbuffzE(:,:,:,2)

    call transpose_y_to_z(this%uhat,zbuff1, this%sp_gpC)
    call this%filzC%filter3(zbuff1,zbuff2,this%nxZ, this%nyZ)
    call transpose_z_to_y(zbuff1,this%uhat, this%sp_gpC)

    call transpose_y_to_z(this%vhat,zbuff1, this%sp_gpC)
    call this%filzC%filter3(zbuff1,zbuff2,this%nxZ, this%nyZ)
    call transpose_z_to_y(zbuff1,this%vhat, this%sp_gpC)

    call transpose_y_to_z(this%what,zbuff3, this%sp_gpE)
    call this%filzE%filter3(zbuff3,zbuff4,this%nxZ, this%nyZ)
    call transpose_z_to_y(zbuff4,this%what, this%sp_gpE)

    if (this%isStratified .or. this%initspinup) then
        call transpose_y_to_z(this%That,zbuff1, this%sp_gpC)
        call this%filzC%filter3(zbuff1,zbuff2,this%nxZ, this%nyZ)
        call transpose_z_to_y(zbuff1,this%That, this%sp_gpC)
    end if

    nullify(zbuff1, zbuff2, zbuff3, zbuff4)
end subroutine

subroutine compute_Sijmean(this, Stmp)
    class(igrid), intent(inout), target :: this
    real(rkind),    dimension(this%gpC%xsz(1),   this%gpC%xsz(2),   this%gpC%xsz(3),   6), intent(out), target :: Stmp

    complex(rkind), dimension(this%sp_gpC%ysz(1),this%sp_gpC%ysz(2),this%sp_gpC%ysz(3),3), target :: SfCtmp
    complex(rkind), dimension(:,:,:), pointer :: u_mean3Dhat, v_mean3Dhat, w_mean3Dhat, ctmpz1, ctmpz2, ctmpy1
    real(rkind),    dimension(:,:,:), pointer :: S11_mean3D, S12_mean3D, S13_mean3D, S22_mean3D, S23_mean3D, S33_mean3D, rbuff1
    real(rkind) :: tidSUMreal

    u_mean3Dhat => SfCtmp(:,:,:,1); v_mean3Dhat => SfCtmp(:,:,:,2); w_mean3Dhat => SfCtmp(:,:,:,3)
    S11_mean3D  => Stmp(:,:,:,1);   S12_mean3D  => Stmp(:,:,:,2);   S13_mean3D  => Stmp(:,:,:,3) 
                                    S22_mean3D  => Stmp(:,:,:,4);   S23_mean3D  => Stmp(:,:,:,5) 
                                                                    S33_mean3D  => Stmp(:,:,:,6) 
    ctmpy1 => this%cbuffyC(:,:,:,1); ctmpz1  => this%cbuffzC(:,:,:,1)
    ctmpz2 => this%cbuffzE(:,:,:,1); rbuff1 => this%rbuffxC(:,:,:,1);


    tidSUMreal = real(this%tidSUM, rkind)

    ! compute forward transforms of mean velocities
    call this%spectC%fft(this%u_mean3D/tidSUMreal, u_mean3Dhat)
    call this%spectC%fft(this%v_mean3D/tidSUMreal, v_mean3Dhat)
    call this%spectC%fft(this%w_mean3D/tidSUMreal, w_mean3Dhat)

    ! dudx
    call this%spectC%mTimes_ik1_oop(u_mean3Dhat, ctmpy1)
    call this%spectC%ifft(ctmpy1, S11_mean3D)
     
    ! dudy
    call this%spectC%mTimes_ik2_oop(u_mean3Dhat, ctmpy1)
    call this%spectC%ifft(ctmpy1, S12_mean3D)
     
    ! dvdx
    call this%spectC%mTimes_ik1_oop(v_mean3Dhat, ctmpy1)
    call this%spectC%ifft(ctmpy1, rbuff1)
    S12_mean3D = half*(S12_mean3D + rbuff1)
     
    ! dvdy
    call this%spectC%mTimes_ik2_oop(v_mean3Dhat, ctmpy1)
    call this%spectC%ifft(ctmpy1, S22_mean3D)
     
    ! dwdx
    call this%spectC%mTimes_ik1_oop(w_mean3Dhat, ctmpy1)
    call this%spectC%ifft(ctmpy1, S13_mean3D)
     
    ! dwdy
    call this%spectC%mTimes_ik2_oop(w_mean3Dhat, ctmpy1)
    call this%spectC%ifft(ctmpy1, S23_mean3D)
    
    ! dudz 
    call transpose_y_to_z(u_mean3Dhat, ctmpz1, this%sp_gpC)
    call this%Pade6opZ%ddz_C2E(ctmpz1, ctmpz2, uBC_bottom, uBC_top)
    call this%Pade6opZ%interpz_E2C(ctmpz2, ctmpz1, dUdzBC_bottom, dUdzBC_top)
    call transpose_z_to_y(ctmpz1, u_mean3Dhat, this%sp_gpC)
    call this%spectC%ifft(u_mean3Dhat, rbuff1)
    S13_mean3D = half*(S13_mean3D + rbuff1)

    ! dvdz 
    call transpose_y_to_z(v_mean3Dhat, ctmpz1, this%sp_gpC)
    call this%Pade6opZ%ddz_C2E(ctmpz1, ctmpz2, vBC_bottom, vBC_top)
    call this%Pade6opZ%interpz_E2C(ctmpz2, ctmpz1, dVdzBC_bottom, dVdzBC_top)
    call transpose_z_to_y(ctmpz1, v_mean3Dhat, this%sp_gpC)
    call this%spectC%ifft(v_mean3Dhat, rbuff1)
    S23_mean3D = half*(S23_mean3D + rbuff1)

    ! dwdz 
    call transpose_y_to_z(w_mean3Dhat, ctmpz1, this%sp_gpC)
    call this%Pade6opZ%ddz_C2E(ctmpz1, ctmpz2, wBC_bottom, wBC_top)
    call this%Pade6opZ%interpz_E2C(ctmpz2, ctmpz1, dWdzBC_bottom, dWdzBC_top)
    call transpose_z_to_y(ctmpz1, w_mean3Dhat, this%sp_gpC)
    call this%spectC%ifft(w_mean3Dhat, S33_mean3D)

    nullify(u_mean3Dhat, v_mean3Dhat, w_mean3Dhat, rbuff1, ctmpy1, ctmpz1, ctmpz2, S11_mean3D, S12_mean3D, S13_mean3D, S22_mean3D, S23_mean3D, S33_mean3D)

end subroutine 


subroutine compute_RapidSlowPressure_Split(this)
     class(igrid), intent(inout), target :: this
     real(rkind), dimension(:,:,:), pointer :: dudx , dudy , dudz , dvdx , dvdy , dvdz , dwdx , dwdy , dwdz

     dudx => this%duidxjC(:,:,:,1); dudy => this%duidxjC(:,:,:,2); dudz => this%duidxjC(:,:,:,3) 
     dvdx => this%duidxjC(:,:,:,4); dvdy => this%duidxjC(:,:,:,5); dvdz => this%duidxjC(:,:,:,6) 
     dwdx => this%duidxjC(:,:,:,7); dwdy => this%duidxjC(:,:,:,8); dwdz => this%duidxjC(:,:,:,9) 
     
     this%prapid = -2.d0*(this%duMdx*dudx + this%duMdy*dvdx + this%duMdz*dwdx &
                     &  + this%dvMdx*dudy + this%dvMdy*dvdy + this%dvMdz*dwdy & 
                     &  + this%dwMdx*dudz + this%dwMdy*dvdz + this%dwMdz*dwdz )  
    
     call this%spectC%fft(this%prapid,this%cbuffyC(:,:,:,1))
     call this%spectC%dealias(this%cbuffyC(:,:,:,1))
     call this%spectC%ifft(this%cbuffyC(:,:,:,1),this%prapid)
     call this%poiss_periodic%poisson_solve(this%prapid)
     this%pslow = this%pressure_dns - this%prapid
     
end subroutine 


subroutine compute_vorticity(this)
     class(igrid), intent(inout), target :: this
     real(rkind),    dimension(:,:,:), pointer :: dudx , dudy , dudzC
     real(rkind),    dimension(:,:,:), pointer :: dvdx , dvdy , dvdzC
     real(rkind),    dimension(:,:,:), pointer :: dwdxC, dwdyC, dwdz
    
     
     dudx  => this%duidxjC(:,:,:,1); dudy  => this%duidxjC(:,:,:,2); dudzC => this%duidxjC(:,:,:,3) 
     dvdx  => this%duidxjC(:,:,:,4); dvdy  => this%duidxjC(:,:,:,5); dvdzC => this%duidxjC(:,:,:,6) 
     dwdxC => this%duidxjC(:,:,:,7); dwdyC => this%duidxjC(:,:,:,8); dwdz  => this%duidxjC(:,:,:,9) 


     this%ox = dwdyC - dvdzC
     this%oy = dudzC - dwdxC
     this%oz = dvdx  - dudy

end subroutine 

subroutine compute_duidxj(this)
    class(igrid), intent(inout), target :: this
    complex(rkind), dimension(:,:,:), pointer :: ctmpz1, ctmpz2
    complex(rkind), dimension(:,:,:), pointer :: ctmpz3, ctmpz4
    complex(rkind), dimension(:,:,:), pointer :: ctmpy1, ctmpy2
    real(rkind),    dimension(:,:,:), pointer :: dudx, dudy, dudz
    real(rkind),    dimension(:,:,:), pointer :: dvdx, dvdy, dvdz
    real(rkind),    dimension(:,:,:), pointer :: dwdx, dwdy, dwdz
    real(rkind),    dimension(:,:,:), pointer :: dvdzC, dudzC 
    real(rkind),    dimension(:,:,:), pointer :: dwdxC, dwdyC
    real(rkind),    dimension(:,:,:), pointer :: dwdzE, dudxE
    real(rkind),    dimension(:,:,:), pointer :: dudyE, dvdxE 
    real(rkind),    dimension(:,:,:), pointer :: dvdyE
    
    complex(rkind), dimension(:,:,:), pointer :: dudxH, dudyH, dudzH 
    complex(rkind), dimension(:,:,:), pointer :: dvdxH, dvdyH, dvdzH
    complex(rkind), dimension(:,:,:), pointer :: dwdxH, dwdyH, dwdzH
    
    complex(rkind), dimension(:,:,:), pointer :: dudxEH, dudyEH, dudzEH 
    complex(rkind), dimension(:,:,:), pointer :: dvdxEH, dvdyEH, dvdzEH
    complex(rkind), dimension(:,:,:), pointer :: dwdxEH, dwdyEH, dwdzEH

    dudx  => this%duidxjC(:,:,:,1); dudy  => this%duidxjC(:,:,:,2); dudzC => this%duidxjC(:,:,:,3) 
    dvdx  => this%duidxjC(:,:,:,4); dvdy  => this%duidxjC(:,:,:,5); dvdzC => this%duidxjC(:,:,:,6) 
    dwdxC => this%duidxjC(:,:,:,7); dwdyC => this%duidxjC(:,:,:,8); dwdz  => this%duidxjC(:,:,:,9) 

    dudxH => this%duidxjChat(:,:,:,1); dudyH => this%duidxjChat(:,:,:,2); dudzH => this%duidxjChat(:,:,:,3) 
    dvdxH => this%duidxjChat(:,:,:,4); dvdyH => this%duidxjChat(:,:,:,5); dvdzH => this%duidxjChat(:,:,:,6) 
    dwdxH => this%duidxjChat(:,:,:,7); dwdyH => this%duidxjChat(:,:,:,8); dwdzH => this%duidxjChat(:,:,:,9) 
    
    dudxEH => this%duidxjEhat(:,:,:,1); dudyEH => this%duidxjEhat(:,:,:,2); dudzEH => this%duidxjEhat(:,:,:,3) 
    dvdxEH => this%duidxjEhat(:,:,:,4); dvdyEH => this%duidxjEhat(:,:,:,5); dvdzEH => this%duidxjEhat(:,:,:,6) 
    dwdxEH => this%duidxjEhat(:,:,:,7); dwdyEH => this%duidxjEhat(:,:,:,8); dwdzEH => this%duidxjEhat(:,:,:,9)
   
    dudxE => this%duidxjE(:,:,:,1); dudyE => this%duidxjE(:,:,:,2); dudz  => this%duidxjE(:,:,:,3)
    dvdxE => this%duidxjE(:,:,:,4); dvdyE => this%duidxjE(:,:,:,5); dvdz  => this%duidxjE(:,:,:,6)
    dwdx  => this%duidxjE(:,:,:,7); dwdy  => this%duidxjE(:,:,:,8); dwdzE => this%duidxjE(:,:,:,9)

    ctmpz1 => this%cbuffzC(:,:,:,1); ctmpz2 => this%cbuffzE(:,:,:,1); 
    ctmpz3 => this%cbuffzC(:,:,:,2); ctmpz4 => this%cbuffzE(:,:,:,2)
    ctmpy1 => this%cbuffyC(:,:,:,1); ctmpy2 => this%cbuffyE(:,:,:,1)

  
    ! dudx
    call this%spectC%mTimes_ik1_oop(this%uhat,dudxH)
    call this%spectC%ifft(dudxH,dudx)
    call this%spectE%mTimes_ik1_oop(this%uEhat,dudxEH)
    call this%spectE%ifft(dudxEH,dudxE)
     
    ! dudy
    call this%spectC%mTimes_ik2_oop(this%uhat,dudyH)
    call this%spectC%ifft(dudyH,dudy)
    call this%spectE%mTimes_ik2_oop(this%uEhat,dudyEH)
    call this%spectE%ifft(dudyEH,dudyE)

    ! dvdx 
    call this%spectC%mTimes_ik1_oop(this%vhat,dvdxH)
    call this%spectC%ifft(dvdxH,dvdx)
    call this%spectE%mTimes_ik1_oop(this%vEhat,dvdxEH)
    call this%spectE%ifft(dvdxEH,dvdxE)

    ! dvdy
    call this%spectC%mTimes_ik2_oop(this%vhat,dvdyH)
    call this%spectC%ifft(dvdyH,dvdy)
    call this%spectE%mTimes_ik2_oop(this%vEhat,dvdyEH)
    call this%spectE%ifft(dvdyEH,dvdyE)

    ! dwdx
    call this%spectC%mTimes_ik1_oop(this%whatC,dwdxH)
    call this%spectC%ifft(dwdxH,dwdxC)
    call this%spectE%mTimes_ik1_oop(this%what, dwdxEH)
    call this%spectE%ifft(dwdxEH,dwdx)

    ! dwdy
    call this%spectC%mTimes_ik2_oop(this%whatC,dwdyH)
    call this%spectC%ifft(dwdyH,dwdyC)
    call this%spectE%mTimes_ik2_oop(this%what,dwdyEH)
    call this%spectE%ifft(dwdyEH,dwdy)
   
    ! dwdz
    call transpose_y_to_z(this%what,ctmpz2,this%sp_gpE)
    call this%Pade6opZ%ddz_E2C(ctmpz2,ctmpz1,wBC_bottom,wBC_top)
    call transpose_z_to_y(ctmpz1,dwdzH,this%sp_gpC)
    call this%spectC%ifft(dwdzH,dwdz)
    call this%Pade6opZ%interpz_C2E(ctmpz1,ctmpz4,dWdzBC_bottom, dWdzBC_top)
    call transpose_z_to_y(ctmpz4,dwdzEH,this%sp_gpE)
    call this%spectE%ifft(dwdzEH,dwdzE)

    ! d2wdz2
    if(.not. this%isinviscid) then
       call this%Pade6opZ%d2dz2_E2E(ctmpz2,ctmpz4,wBC_bottom,wBC_top)
       call transpose_z_to_y(ctmpz4,this%d2wdz2hatE,this%sp_gpE)
    end if

    ! dudz and d2udz2
    if ((topWall == 1).or.(botWall == 1)) then
        call transpose_y_to_z(this%uEhat,ctmpz2,this%sp_gpE)
        ! Compute dudzC
        call this%Pade6opZ%ddz_E2C(ctmpz2,ctmpz1,uBC_bottom,uBC_top)
        call transpose_z_to_y(ctmpz1,dudzH,this%sp_gpC)
        call this%spectC%ifft(dudzH,dudzC)
        ! Now compute d2udz2C
        call this%Pade6opZ%ddz_C2C(ctmpz1,ctmpz3,-uBC_bottom,-uBC_top)
        call transpose_z_to_y(ctmpz3,this%d2udz2hatC,this%sp_gpC)
        ! Now compute dudzE
        call this%Pade6opZ%interpz_C2E(ctmpz1,ctmpz4,-uBC_bottom,-uBC_top) 
        call transpose_z_to_y(ctmpz4,dudzEH,this%sp_gpC)
        call this%spectE%ifft(dudzEH,dudz)
    else
        call transpose_y_to_z(this%uhat,ctmpz1,this%sp_gpC)
        call this%Pade6opZ%ddz_C2E(ctmpz1,ctmpz2,uBC_bottom,uBC_top)
        call transpose_z_to_y(ctmpz2,dudzEH,this%sp_gpE)
        call this%spectE%ifft(dudzEH,dudz)
        if (.not. this%isinviscid) then
               if ((uBC_top == 0) .or. (uBC_bottom == 0)) then
                  call this%Pade6opZ%ddz_C2E(ctmpz1,ctmpz4,uBC_bottom,uBC_top)
                  call this%Pade6opZ%ddz_E2C(ctmpz4,ctmpz3,dUdzBC_bottom,dUdzBC_top)
               else
                  call this%Pade6opZ%d2dz2_C2C(ctmpz1,ctmpz3,uBC_bottom,uBC_top)
               end if
               call transpose_z_to_y(ctmpz3,this%d2udz2hatC,this%sp_gpC)
        end if 
        call this%Pade6opZ%interpz_E2C(ctmpz2,ctmpz1,dUdzBC_bottom,dUdzBC_top)
        call transpose_z_to_y(ctmpz1,dudzH,this%sp_gpC)
        call this%spectC%ifft(dudzH,dudzC)
    end if 

    ! dvdz and d2vdz2
    if ((topWall == 1).or.(botWall == 1)) then
        call transpose_y_to_z(this%vEhat,ctmpz2,this%sp_gpE)
        ! compute dvdzC
        call this%Pade6opZ%ddz_E2C(ctmpz2,ctmpz1,vBC_bottom,vBC_top)
        call transpose_z_to_y(ctmpz1,dvdzH,this%sp_gpC)
        call this%spectC%ifft(dvdzH,dvdzC)
        ! Now compute d2udz2C
        call this%Pade6opZ%ddz_C2C(ctmpz1,ctmpz3,-vBC_bottom,-vBC_top)
        call transpose_z_to_y(ctmpz3,this%d2vdz2hatC,this%sp_gpC)
        ! Now compute dvdzE
        call this%Pade6opZ%interpz_C2E(ctmpz1,ctmpz4,-vBC_bottom,-vBC_top) 
        call transpose_z_to_y(ctmpz4,dvdzEH,this%sp_gpC)
        call this%spectE%ifft(dvdzEH,dvdz)
    else
        call transpose_y_to_z(this%vhat,ctmpz1,this%sp_gpC)
        call this%Pade6opZ%ddz_C2E(ctmpz1,ctmpz2,vBC_bottom,vBC_top)
        call transpose_z_to_y(ctmpz2,dvdzEH,this%sp_gpE)
        call this%spectE%ifft(dvdzEH,dvdz)
        if (.not. this%isinviscid) then
            if ((vBC_top == 0) .or. (vBC_bottom == 0)) then
               call this%Pade6opZ%ddz_C2E(ctmpz1,ctmpz4,vBC_bottom,vBC_top)
               call this%Pade6opZ%ddz_E2C(ctmpz4,ctmpz3,dVdzBC_bottom,dVdzBC_top)
            else
               call this%Pade6opZ%d2dz2_C2C(ctmpz1,ctmpz3,vBC_bottom,vBC_top)
            end if
            call transpose_z_to_y(ctmpz3,this%d2vdz2hatC,this%sp_gpC)
        end if 
        call this%Pade6opZ%interpz_E2C(ctmpz2,ctmpz1,dVdzBC_bottom,dVdzBC_top)
        call transpose_z_to_y(ctmpz1,dvdzH,this%sp_gpC)
        call this%spectC%ifft(dvdzH,dvdzC)
    end if 

end subroutine


subroutine compute_dTdxi(this)
    class(igrid), intent(inout), target :: this
    complex(rkind), dimension(:,:,:), pointer :: ctmpz1, ctmpz2, ctmpz3
    complex(rkind), dimension(:,:,:), pointer :: ctmpy1

    ctmpz1 => this%cbuffzC(:,:,:,1); ctmpz2 => this%cbuffzE(:,:,:,1); 
    ctmpy1 => this%cbuffyC(:,:,:,1); ctmpz3 => this%cbuffzC(:,:,:,2)

    call this%spectC%mtimes_ik1_oop(this%That,this%dTdxH)
    call this%spectC%ifft(this%dTdxH,this%dTdxC)

    call this%spectC%mtimes_ik2_oop(this%That,this%dTdyH)
    call this%spectC%ifft(this%dTdyH,this%dTdyC)

    call this%spectE%mtimes_ik1_oop(this%TEhat,this%cbuffyE(:,:,:,1))
    call this%spectE%ifft(this%cbuffyE(:,:,:,1),this%dTdxE)
    
    call this%spectE%mtimes_ik2_oop(this%TEhat,this%cbuffyE(:,:,:,1))
    call this%spectE%ifft(this%cbuffyE(:,:,:,1),this%dTdyE)
   
    if (((this%botBC_Temp == 3) .and. (this%topBC_Temp == 3)).or.(this%botBC_Temp == 2)) then ! symmetric dirichlet BCs
        call transpose_y_to_z(this%TEhat, ctmpz2, this%sp_gpE)
        call this%Pade6opZ%ddz_E2C(ctmpz2,ctmpz1,TBC_bottom,TBC_top)
        call transpose_z_to_y(ctmpz1,this%dTdzHC,this%sp_gpC)
        call this%spectC%ifft(this%dTdzHC,this%dTdzC)
        
        call this%Pade6opZ%interpz_C2E(ctmpz1,ctmpz2,dTdzBC_bottom,dTdzBC_top)
        call transpose_z_to_y(ctmpz2, this%dTdzH, this%sp_gpE)
        call this%spectE%ifft(this%dTdzH,this%dTdzE)

        if (.not. this%isInviscid) then
            call this%Pade6opZ%ddz_C2C(ctmpz1,ctmpz3,dTdzBC_bottom,dTdzBC_top)  
            call transpose_z_to_y(ctmpz3,this%d2Tdz2hatC,this%sp_gpC)
        end if 
        
    else
        call transpose_y_to_z(this%That, ctmpz1, this%sp_gpC)
        call this%Pade6opZ%ddz_C2E(ctmpz1,ctmpz2,TBC_bottom,TBC_top)
        if (.not. this%isInviscid) then
            call this%Pade6opZ%d2dz2_C2C(ctmpz1,ctmpz3,TBC_bottom, TBC_top)    
            call transpose_z_to_y(ctmpz3,this%d2Tdz2hatC,this%sp_gpC)
        end if 
        call this%Pade6opZ%interpz_E2C(ctmpz2,ctmpz1,dTdzBC_bottom,dTdzBC_top)

        call transpose_z_to_y(ctmpz2, this%dTdzH, this%sp_gpE)
        call this%spectE%ifft(this%dTdzH,this%dTdzE)
        
        call transpose_z_to_y(ctmpz1,this%dTdzHC,this%sp_gpC)
        call this%spectC%ifft(this%dTdzHC,this%dTdzC)
    end if 

end subroutine
