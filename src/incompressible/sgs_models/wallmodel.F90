subroutine destroyWallModel(this)
   class(sgs_igrid), intent(inout) :: this
   deallocate(this%tauijWM, this%tauijWMhat_inZ, this%tauijWMhat_inY)
   if (allocated(this%filteredSpeedSq)) deallocate(this%filteredSpeedSq)
end subroutine

subroutine initWallModel(this, SurfaceFilterFact)
   class(sgs_igrid), intent(inout) :: this
   real(rkind), intent(in) :: SurfaceFilterFact

   this%useWallModel = .true.
   allocate(this%tauijWM(this%gpE%xsz(1),this%gpE%xsz(2),this%gpE%xsz(3),2))
   allocate(this%tauijWMhat_inZ(this%sp_gpE%zsz(1),this%sp_gpE%zsz(2),this%sp_gpE%zsz(3),2))
   allocate(this%tauijWMhat_inY(this%sp_gpE%ysz(1),this%sp_gpE%ysz(2),this%sp_gpE%ysz(3),2))
   this%tauijWM = 0.d0
   this%tauijWMhat_inZ = dcmplx(0.d0, 0.d0)
   this%tauijWMhat_inY = dcmplx(0.d0, 0.d0)
   select case(this%WallModel)
   case (1) ! Standard Moeng Wall model

   case (2) ! Bou-Zeid Wall model
      allocate(this%filteredSpeedSq(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)))
   case default
      call gracefulExit("Invalid choice of Wallmodel.",324)
   end select

   if (this%useFullyLocalWM) then
      allocate(this%usurf_filt(this%gpC%zsz(1),this%gpC%zsz(2)))
      allocate(this%vsurf_filt(this%gpC%zsz(1),this%gpC%zsz(2)))
      allocate(this%Tmatch_filt(this%gpC%zsz(1),this%gpC%zsz(2)))
      allocate(this%ustar_surf(this%gpC%zsz(1),this%gpC%zsz(2)))
      allocate(this%wTheta_surf(this%gpC%zsz(1),this%gpC%zsz(2)))
      allocate(this%PsiM_surf(this%gpC%zsz(1),this%gpC%zsz(2)))
      allocate(this%Linv_surf(this%gpC%zsz(1),this%gpC%zsz(2)))
      allocate(this%T_surf(this%gpC%zsz(1),this%gpC%zsz(2)))
      allocate(this%filteredSpeedSq(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3))) ! Howland: Added 1/25/21
      call this%spectC%ResetSurfaceFilter(SurfaceFilterFact)
      call message(2,"Fully local wall model set up with a filter factor:", SurfaceFilterFact)
   end if 

   if (this%isStratified .or. this%initspinup) then
      allocate(this%q3HAT_AtWall(this%sp_gpE%zsz(1),this%sp_gpE%zsz(2)))
      this%q3HAT_AtWall = dcmplx(0.d0, 0.d0)
   end if 
   
   this%T_surf_mean=0.d0
   
end subroutine


subroutine computeWallStress(this, u, v, T, uhat, vhat, That)
   class(sgs_igrid), intent(inout) :: this
   complex(rkind), dimension(this%sp_gpC%ysz(1),this%sp_gpC%ysz(2),this%sp_gpC%ysz(3)), intent(in) :: uhat, vhat, That
   real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)), intent(in) :: u, v, T
    complex(rkind), dimension(:,:,:), pointer :: cbuffz, cbuffy
  
    
   cbuffz => this%cbuffzC(:,:,:,1)
   cbuffy => this%cbuffyC(:,:,:,1)
   
   call this%compute_and_bcast_surface_Mn(u, v, uhat, vhat, That)
   
   if (this%useFullyLocalWM) then
        ! No temperature filter
        !call this%getfilteredMatchingVelocity(uhat, vhat, T)
        ! Filter temperature
        call this%getfilteredMatchingVelocity(uhat, vhat, That)
        call this%compute_surface_stress()
   else
        
        select case (this%WallModel)
        case (1) ! Standard Moeng Wall model
           this%WallMFactor = -this%ustar*this%ustar/(this%Uspmn + 1.D-13)

           ! Tau_13
           call transpose_y_to_z(uhat, cbuffz, this%sp_gpC)
           this%tauijWMhat_inZ(:,:,1,1) = this%WallMFactor*cbuffz(:,:,this%WM_matchingIndex) 
           call transpose_z_to_y(this%tauijWMhat_inZ(:,:,:,1), this%tauijWMhat_inY(:,:,:,1), this%sp_gpE)
           call this%spectE%ifft(this%tauijWMhat_inY(:,:,:,1), this%tauijWM(:,:,:,1))

           ! Tau_23
           call transpose_y_to_z(vhat, cbuffz, this%sp_gpC)
           this%tauijWMhat_inZ(:,:,1,2) = this%WallMFactor*cbuffz(:,:,this%WM_matchingIndex) 
           call transpose_z_to_y(this%tauijWMhat_inZ(:,:,:,2), this%tauijWMhat_inY(:,:,:,2), this%sp_gpE)
           call this%spectE%ifft(this%tauijWMhat_inY(:,:,:,2), this%tauijWM(:,:,:,2))

        case (2) ! Bou-zeid Wall model 
           this%WallMFactor = -(kappa/(log(this%dz/(two*this%z0)) - this%PsiM))**2 
           call this%getfilteredSpeedSqAtWall(uhat, vhat)
           
           call this%spectC%fft(this%filteredSpeedSq, cbuffy)
           call transpose_y_to_z(cbuffy, cbuffz, this%sp_gpC)
           
           ! tau_13
           this%tauijWMhat_inZ(:,:,1,1) = (this%WallMFactor*this%umn/this%Uspmn) * cbuffz(:,:,this%WM_matchingIndex) 
           call transpose_z_to_y(this%tauijWMhat_inZ(:,:,:,1), this%tauijWMhat_inY(:,:,:,1), this%sp_gpE)
           call this%spectE%ifft(this%tauijWMhat_inY(:,:,:,1), this%tauijWM(:,:,:,1))
           
           ! tau_23
           this%tauijWMhat_inZ(:,:,1,2) = (this%WallMFactor*this%vmn/this%Uspmn) * cbuffz(:,:,this%WM_matchingIndex) 
           call transpose_z_to_y(this%tauijWMhat_inZ(:,:,:,2), this%tauijWMhat_inY(:,:,:,2), this%sp_gpE)
           call this%spectE%ifft(this%tauijWMhat_inY(:,:,:,2), this%tauijWM(:,:,:,2))
        end select
   end if 
end subroutine

subroutine embed_WM_stress(this)
   class(sgs_igrid), intent(inout) :: this
   
   if(this%gpE%xst(3)==1) then
      this%tau_13(:,:,1) = this%tauijWM(:,:,1,1)
      this%tau_23(:,:,1) = this%tauijWM(:,:,1,2)
   endif
end subroutine 

subroutine embed_WM_PotTflux(this)
   class(sgs_igrid), intent(inout) :: this

   if(this%gpE%xst(3)==1) then
      this%q3E(:,:,1) = this%wTh_surf
   endif
end subroutine 


subroutine computeWall_PotTFlux(this)
   class(sgs_igrid), intent(inout) :: this
 
   if (nrank == 0) then
      this%q3HAT_AtWall(1,1) = cmplx(this%wTh_surf/this%MeanFact,zero,rkind)
   end if

end subroutine

subroutine compute_and_bcast_surface_Mn(this, u, v, uhat, vhat, That )
    use mpi
    !use constants, only: four
    use kind_parameters, only: mpirkind
    class(sgs_igrid), intent(inout), target :: this
    complex(rkind), dimension(this%sp_gpC%ysz(1),this%sp_gpC%ysz(2),this%sp_gpC%ysz(3)), intent(in) :: uhat, vhat, That
    real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)), intent(in) :: u, v
    real(rkind), dimension(:,:,:), pointer :: rbuff
    complex(rkind), dimension(:,:,:), pointer :: cbuff
    integer :: ierr

    rbuff => this%rbuffxC(:,:,:,1)
    cbuff => this%cbuffyC(:,:,:,1)
    rbuff = sqrt(u*u + v*v)
    call this%spectC%fft(rbuff,cbuff)

    if (nrank == 0) then
        this%Umn = real(uhat(1,1,this%WM_matchingIndex),rkind)*this%meanFact
        this%Vmn = real(vhat(1,1,this%WM_matchingIndex),rkind)*this%meanFact
        this%Uspmn = real(cbuff(1,1,this%WM_matchingIndex),rkind)*this%meanFact
        if (this%isStratified .or. this%initSpinUp) this%Tmn = real(That(1,1,1),rkind)*this%meanFact
    end if
    call mpi_bcast(this%Umn,1,mpirkind,0,mpi_comm_world,ierr)
    call mpi_bcast(this%Vmn,1,mpirkind,0,mpi_comm_world,ierr)
    call mpi_bcast(this%Uspmn,1,mpirkind,0,mpi_comm_world,ierr)
    if (this%isStratified .or. this%initSpinup) call mpi_bcast(this%Tmn,1,mpirkind,0,mpi_comm_world,ierr)

    call this%getSurfaceQuantities() 
end subroutine

subroutine getfilteredSpeedSqAtWall(this, uhatC, vhatC)
    class(sgs_igrid), intent(inout), target :: this
    complex(rkind), dimension(this%sp_gpC%ysz(1),this%sp_gpC%ysz(2),this%sp_gpC%ysz(3)), intent(in) :: uhatC, vhatC

    real(rkind), dimension(:,:,:), pointer :: rbuffx1, rbuffx2
    complex(rkind), dimension(:,:,:), pointer :: cbuffy, tauWallH

    cbuffy => this%cbuffyC(:,:,:,1); tauWallH => this%cbuffzC(:,:,:,1)     
    rbuffx1 => this%filteredSpeedSq; rbuffx2 => this%rbuffxC(:,:,:,1)

    call transpose_y_to_z(uhatC,tauWallH,this%sp_gpC)
    call this%spectC%SurfaceFilter_ip(tauWallH(:,:,1))
    call transpose_z_to_y(tauWallH,cbuffy, this%sp_gpC)
    call this%spectC%ifft(cbuffy,rbuffx1)

    call transpose_y_to_z(vhatC,tauWallH,this%sp_gpC)
    call this%spectC%SurfaceFilter_ip(tauWallH(:,:,1))
    call transpose_z_to_y(tauWallH,cbuffy, this%sp_gpC)
    call this%spectC%ifft(cbuffy,rbuffx2)

    rbuffx1 = rbuffx1*rbuffx1
    rbuffx2 = rbuffx2*rbuffx2
    rbuffx1 = rbuffx1 + rbuffx2

end subroutine  

subroutine getfilteredMatchingVelocity(this, uhatC, vhatC, That)
    class(sgs_igrid), intent(inout), target :: this
    complex(rkind), dimension(this%sp_gpC%ysz(1),this%sp_gpC%ysz(2),this%sp_gpC%ysz(3)), intent(in) :: uhatC, vhatC, That
    !real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)), intent(in) :: T
    real(rkind), dimension(:,:,:), pointer :: rbuffx1, rbuffx2
    complex(rkind), dimension(:,:,:), pointer :: cbuffy, tauWallH

    cbuffy => this%cbuffyC(:,:,:,1); tauWallH => this%cbuffzC(:,:,:,1)     
    rbuffx1 => this%filteredSpeedSq; rbuffx2 => this%rbuffxC(:,:,:,1)
    
    call transpose_y_to_z(uhatC,tauWallH,this%sp_gpC)
    call this%spectC%SurfaceFilter_ip(tauWallH(:,:,1))
    call transpose_z_to_y(tauWallH,cbuffy, this%sp_gpC)
    call this%spectC%ifft(cbuffy,rbuffx1)
    call transpose_x_to_y(rbuffx1,this%rbuffyC(:,:,:,1),this%gpC)
    call transpose_y_to_z(this%rbuffyC(:,:,:,1),this%rbuffzC(:,:,:,1), this%gpC)
    this%usurf_filt = this%rbuffzC(:,:,1,1) 

    call transpose_y_to_z(vhatC,tauWallH,this%sp_gpC)
    call this%spectC%SurfaceFilter_ip(tauWallH(:,:,1))
    call transpose_z_to_y(tauWallH,cbuffy, this%sp_gpC)
    call this%spectC%ifft(cbuffy,rbuffx2)
    call transpose_x_to_y(rbuffx2,this%rbuffyC(:,:,:,1),this%gpC)
    call transpose_y_to_z(this%rbuffyC(:,:,:,1),this%rbuffzC(:,:,:,1), this%gpC)
    this%vsurf_filt = this%rbuffzC(:,:,1,1) 

    if (this%isStratified) then
        ! Filter for temperature
        call transpose_y_to_z(That,tauWallH,this%sp_gpC)
        call this%spectC%SurfaceFilter_ip(tauWallH(:,:,1))
        call transpose_z_to_y(tauWallH,cbuffy, this%sp_gpC)
        call this%spectC%ifft(cbuffy,rbuffx1)
        call transpose_x_to_y(rbuffx1,this%rbuffyC(:,:,:,1),this%gpC)
        call transpose_y_to_z(this%rbuffyC(:,:,:,1),this%rbuffzC(:,:,:,1), this%gpC)
        this%Tmatch_filt = this%rbuffzC(:,:,1,1) 
        
        ! No filter for temperature
        !call transpose_x_to_y(T,this%rbuffyC(:,:,:,1),this%gpC)
        !call transpose_y_to_z(this%rbuffyC(:,:,:,1),this%rbuffzC(:,:,:,1), this%gpC)
        !this%Tmatch_filt = this%rbuffzC(:,:,1,1) 
    else
        this%Tmatch_filt = 0.d0 
    end if 
end subroutine 

subroutine compute_surface_stress(this)
    class(sgs_igrid), intent(inout) :: this
    integer :: i, j

    ! Compute local wall-quantities at each wall node
    do j = 1,this%gpC%zsz(2)
        do i = 1,this%gpC%zsz(1)
            call this%compute_local_wallmodel(this%usurf_filt(i,j), this%vsurf_filt(i,j), this%Tmatch_filt(i,j), & 
                this%wTheta_surf(i,j), this%ustar_surf(i,j), this%Linv_surf(i,j), this%PsiM_surf(i,j), this%T_surf(i,j))
        end do 
    end do 

    ! Compute stress and heat flux 
    this%rbuffzE(:,:,1,1) = -(this%ustar_surf**2)*(this%usurf_filt)/(sqrt(this%usurf_filt**2 + this%vsurf_filt**2) + 1.d-13)
    call transpose_z_to_y(this%rbuffzE(:,:,:,1), this%rbuffyE(:,:,:,1), this%gpE)
    call transpose_y_to_x(this%rbuffyE(:,:,:,1), this%tauijWM(:,:,:,1), this%gpE)
    
    this%rbuffzE(:,:,1,1) = -(this%ustar_surf**2)*(this%vsurf_filt)/(sqrt(this%usurf_filt**2 + this%vsurf_filt**2) + 1.d-13)
    call transpose_z_to_y(this%rbuffzE(:,:,:,1), this%rbuffyE(:,:,:,1), this%gpE)
    call transpose_y_to_x(this%rbuffyE(:,:,:,1), this%tauijWM(:,:,:,2), this%gpE)

    ! Surface heat flux (horizontally averaged can be avoided) 
    this%wTh_surf = p_sum(sum(this%wTheta_surf))/real(this%gpC%xsz(1)*this%gpC%ysz(2),rkind)
    this%Tsurf = p_sum(sum(this%T_surf))/real(this%gpC%xsz(1)*this%gpC%ysz(2),rkind)

    ! Store horizontally averaged things
    this%ustar       = p_sum(sum(this%ustar_surf))/real(this%gpC%xsz(1)*this%gpC%ysz(2),rkind)
    this%InvObLength = p_sum(sum(this%Linv_surf))/real(this%gpC%xsz(1)*this%gpC%ysz(2),rkind)
end subroutine 

subroutine compute_local_wallmodel(this, ux, uy, Tmn, wTh_surf, ustar, Linv, PsiM, T_surf)
    class(sgs_igrid), intent(inout) :: this
    real(rkind), intent(in) :: ux, uy, Tmn
    real(rkind), intent(out) :: wTh_surf, ustar, Linv, PsiM, T_surf

    real(rkind) :: ustarNew, ustarDiff, dTheta, at
    real(rkind) :: a, b, c, PsiH, wTh, u, xi, xisq
    real(rkind) :: hwm
    integer, parameter :: itermax = 100 
    integer :: idx

    hwm = this%dz/two 
    if (this%isStratified) then
      select case (this%botBC_Temp)
      case(0) ! Dirichlet BC for temperature 
          dTheta = this%Tsurf - Tmn; Linv = zero
          ustarDiff = one; wTh = zero
          a=log(hwm/this%z0); b=beta_h*hwm; c=beta_m*hwm
          PsiM = zero; PsiH = zero; idx = 0; ustar = one; u = sqrt(ux*ux + uy*uy)
          at=log(hwm/this%z0t)

          do while ( (ustarDiff > 1d-12) .and. (idx < itermax))
              ustarNew = u*kappa/(a - PsiM)
              wTh = dTheta*ustarNew*kappa/(at - PsiH) 
              Linv = -kappa*wTh/((this%Fr**2) * this%ThetaRef*ustarNew**3)
              if (this%WallFunctionType == 2) then
                call this%getMO_wallfunction(hwm, Linv, PsiM, PsiH)
              else
                if (Linv .ge. zero) then 
                  ! similarity functions if stable stratification is present
                  PsiM = -c*Linv;         PsiH = -b*Linv; 
                else
                  ! similarity functions if unstable stratification is present
                  xisq = sqrt(one-15.d0*hwm*Linv); xi = sqrt(xisq)
                  PsiM = two*log(half*(one+xi)) + log(half*(one+xisq)) - two*atan(xi) + piby2; 
                  PsiH = two*log(half*(one+xisq));
                endif
              end if 
              ustarDiff = abs((ustarNew - ustar)/ustarNew)
              ustar = ustarNew; idx = idx + 1
          end do 
          wTh_surf = wTh
          T_surf = this%Tsurf
      case(1) ! Homogeneous Neumann BC for temperature
          ustar = this%Uspmn*kappa/(log(hwm/this%z0))
          Linv = zero
          wTh_surf = zero
          PsiM = zero
          T_surf = Tmn 
      case(2) ! Inhomogeneous Neumann BC for temperature
          Linv = zero; !dTheta = this%Tsurf - this%Tmn;
          ustarDiff = one; wTh = this%wTh_surf
          a=log(hwm/this%z0); b=beta_h*hwm; c=beta_m*hwm
          PsiM = zero; PsiH = zero; idx = 0; ustar = one; u = sqrt(ux*ux + uy*uy) 
          at=log(hwm/this%z0t)
   
          do while ( (ustarDiff > 1d-12) .and. (idx < itermax))
              ustarNew = u*kappa/(a - PsiM)
              Linv = -kappa*wTh/((this%Fr**2) * this%ThetaRef*ustarNew**3)
              if (this%WallFunctionType == 2) then
                call this%getMO_wallfunction(hwm, Linv, PsiM, PsiH)
              else
                if (Linv .ge. zero) then 
                  ! similarity functions if stable stratification is present
                  PsiM = -c*Linv;         PsiH = -b*Linv; 
                else
                  ! similarity functions if unstable stratification is present
                  xisq = sqrt(one-15.d0*hwm*Linv); xi = sqrt(xisq)
                  PsiM = two*log(half*(one+xi)) + log(half*(one+xisq)) - two*atan(xi) + piby2; 
                  PsiH = two*log(half*(one+xisq));
                endif
              end if 
              ustarDiff = abs((ustarNew - ustar)/ustarNew)
              ustar = ustarNew; idx = idx + 1
          end do
          wTh_surf = this%wTh_surf
          T_surf = Tmn + wTh*(at-PsiH)/(ustar*kappa)
      end select
   else
          ustar = sqrt(ux*ux + uy*uy)*kappa/(log(hwm/this%z0))
          Linv = zero
          wTh_surf = zero
          PsiM = zero
          T_surf = Tmn 
    end if
 
   this%T_surf_mean = T_surf !p_sum(sum(T_surf))/real(this%gpC%xsz(1)*this%gpC%ysz(2),rkind)
end subroutine 

subroutine getMO_wallfunction(this, z, Linv, PsiM, PsiH)
    class(sgs_igrid), intent(inout) :: this
    real(rkind), intent(in) :: z, Linv
    real(rkind), intent(out) :: PsiM, PsiH
    real(rkind) :: zL
    real(rkind), parameter :: a = 2.5d0, b = 1.1d0  

    zL = z*Linv

    if (Linv < 0.d0) then
        PsiM = (1.d0 - 15.2d0*zL)**(-0.25d0)
        PsiH = (1.d0 - 15.2d0*zL)**(-0.5d0)
    else
        PsiM = 1.d0 + 6.1d0*(zL + (zL**a) * (1.d0 + zL**a)**(-1.d0 + 1.d0/a) )/(zL + (1.d0 + zL**a)**(1.d0/a))
        PsiH = 1.d0 + 5.3d0*(zL + (zL**b) * (1.d0 + zL**b)**(-1.d0 + 1.d0/b) )/(zL + (1.d0 + zL**b)**(1.d0/b))
    end if 

end subroutine 

subroutine getSurfaceQuantities(this)
    class(sgs_igrid), intent(inout) :: this
    integer :: idx
    integer, parameter :: itermax = 100 
    real(rkind) :: ustarNew, ustarDiff, dTheta, ustar, at
    real(rkind) :: a, b, c, PsiH, PsiM, wTh, u, Linv, xi, xisq
    real(rkind) :: hwm

    hwm = this%dz/two + (this%WM_matchingIndex - 1)*this%dz
    if (this%isStratified) then
      select case (this%botBC_Temp)
      case(0) ! Dirichlet BC for temperature 
          dTheta = this%Tsurf - this%Tmn; Linv = zero
          ustarDiff = one; wTh = zero
          a=log(hwm/this%z0); b=beta_h*hwm; c=beta_m*hwm
          PsiM = zero; PsiH = zero; idx = 0; ustar = one; u = this%Uspmn
          at=log(hwm/this%z0t)

          do while ( (ustarDiff > 1d-12) .and. (idx < itermax))
              ustarNew = u*kappa/(a - PsiM)
              wTh = dTheta*ustarNew*kappa/(at - PsiH) 
              Linv = -kappa*wTh/((this%Fr**2) * this%ThetaRef*ustarNew**3)
              if (Linv .ge. zero) then 
                ! similarity functions if stable stratification is present
                PsiM = -c*Linv;         PsiH = -b*Linv; 
              else
                ! similarity functions if unstable stratification is present
                xisq = sqrt(one-15.d0*hwm*Linv); xi = sqrt(xisq)
                PsiM = two*log(half*(one+xi)) + log(half*(one+xisq)) - two*atan(xi) + piby2; 
                PsiH = two*log(half*(one+xisq));
              endif
              ustarDiff = abs((ustarNew - ustar)/ustarNew)
              ustar = ustarNew; idx = idx + 1
          end do 
          this%ustar = ustar; this%invObLength = Linv; this%wTh_surf = wTh
          this%PsiM = PsiM
      case(1) ! Homogeneous Neumann BC for temperature
          this%ustar = this%Uspmn*kappa/(log(hwm/this%z0))
          this%invObLength = zero
          this%wTh_surf = zero
          this%PsiM = zero
      case(2) ! Inhomogeneous Neumann BC for temperature
          Linv = zero; !dTheta = this%Tsurf - this%Tmn;
          ustarDiff = one; wTh = this%wTh_surf
          a=log(hwm/this%z0); b=beta_h*hwm; c=beta_m*hwm
          PsiM = zero; PsiH = zero; idx = 0; ustar = one; u = this%Uspmn
          at=log(hwm/this%z0t)
   
          do while ( (ustarDiff > 1d-12) .and. (idx < itermax))
              ustarNew = u*kappa/(a - PsiM)
              Linv = -kappa*wTh/((this%Fr**2) * this%ThetaRef*ustarNew**3)
              if (Linv .ge. zero) then 
                ! similarity functions if stable stratification is present
                PsiM = -c*Linv;         PsiH = -b*Linv; 
              else
                ! similarity functions if unstable stratification is present
                xisq = sqrt(one-15.d0*hwm*Linv); xi = sqrt(xisq)
                PsiM = two*log(half*(one+xi)) + log(half*(one+xisq)) - two*atan(xi) + piby2; 
                PsiH = two*log(half*(one+xisq));
              endif
              ustarDiff = abs((ustarNew - ustar)/ustarNew)
              ustar = ustarNew; idx = idx + 1
          end do 
          this%ustar = ustar; this%invObLength = Linv; 
          this%Tsurf = this%Tmn + wTh*(at-PsiH)/(ustar*kappa)
          this%PsiM = PsiM
      end select
   else
          this%ustar = this%Uspmn*kappa/(log(hwm/this%z0))
          this%invObLength = zero
          this%wTh_surf = zero
          this%PsiM = zero
    end if 
end subroutine
