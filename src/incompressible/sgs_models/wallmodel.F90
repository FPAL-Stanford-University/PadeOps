subroutine destroyWallModel(this)
   class(sgs_igrid), intent(inout) :: this
   deallocate(this%tauijWM, this%tauijWMhat_inZ, this%tauijWMhat_inY)
   if (allocated(this%filteredSpeedSq)) deallocate(this%filteredSpeedSq)
end subroutine

subroutine initWallModel(this)
   class(sgs_igrid), intent(inout) :: this

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

   if (this%isStratified .or. this%initspinup) then
      allocate(this%q3HAT_AtWall(this%sp_gpE%zsz(1),this%sp_gpE%zsz(2)))
      this%q3HAT_AtWall = dcmplx(0.d0, 0.d0)
   end if 
end subroutine


subroutine computeWallStress(this, u, v, uhat, vhat, That)
   class(sgs_igrid), intent(inout) :: this
   complex(rkind), dimension(this%sp_gpC%ysz(1),this%sp_gpC%ysz(2),this%sp_gpC%ysz(3)), intent(in) :: uhat, vhat, That
   real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)), intent(in) :: u, v
    complex(rkind), dimension(:,:,:), pointer :: cbuffz, cbuffy
  
   cbuffz => this%cbuffzC(:,:,:,1)
   cbuffy => this%cbuffyC(:,:,:,1)
   call this%compute_and_bcast_surface_Mn(u, v, uhat, vhat, That)
   
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
      call this%getfilteredSpeedSqAtWall(uhat, vhat)
      if(this%is_z0_varying) then
          this%WallMFactorvar = -(kappa/(log(this%dz/(two*this%z0var)) - this%PsiM))**2 
          call this%BouZeidLocalModel()
      else 
          this%WallMFactor = -(kappa/(log(this%dz/(two*this%z0)) - this%PsiM))**2 
          
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
      endif

   case (3) ! Abkar-PA (2012) heterogeneous model
      !this%kaplnzfac = this%kappa/(half*this%dz/this%z0s)
      !call this%getSpanAvgVelAtWall(uhat, vhat)
      !where(this%lamfact < (1.0d0-1.0d-5) )
      !    this%ustarsqvar = this%SpanAvgSpeed - this%lamfact*ust1fac
      !    this%ustarsqvar = this%ustarsqvar/(one - this%lamfact) * kaplnzfac
      !else
      !    this%ustarsqvar = this%SpanAvgSpeed * kaplnzfac
      !end where
   end select

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

subroutine BouZeidLocalModel(this)
    class(sgs_igrid), intent(inout), target :: this
    real(rkind), dimension(:,:,:), pointer :: rbuffx1, rbuffx2, rbuffx3

    rbuffx1 => this%rbuffxC(:,:,:,1);    rbuffx2 => this%rbuffxC(:,:,:,2)
    rbuffx3 => this%rbuffxC(:,:,:,3);

    rbuffx3(:,:,1) = sqrt(this%filteredSpeedSq(:,:,1))

    ! tau_13
    rbuffx1(:,:,1) = this%WallmFactorvar * this%Uxvar * rbuffx3(:,:,1)
    call transpose_x_to_y(rbuffx1, this%rbuffyC(:,:,:,1), this%gpC)
    call transpose_y_to_z(this%rbuffyC(:,:,:,1), this%rbuffzC(:,:,:,1), this%gpC)

    this%rbuffzE = 0.0d0;    this%rbuffzE(:,:,1,1) = this%rbuffzC(:,:,1,1)
    call transpose_z_to_y(this%rbuffzE(:,:,:,1), this%rbuffyE(:,:,:,1), this%gpE)
    call transpose_y_to_x(this%rbuffyE(:,:,:,1), this%tauijWM(:,:,:,1), this%gpE)

    ! tau_23
    rbuffx2(:,:,1) = this%WallmFactorvar * this%Uyvar * rbuffx3(:,:,1)
    call transpose_x_to_y(rbuffx2, this%rbuffyC(:,:,:,1), this%gpC)
    call transpose_y_to_z(this%rbuffyC(:,:,:,1), this%rbuffzC(:,:,:,1), this%gpC)

    this%rbuffzE = 0.0d0;    this%rbuffzE(:,:,1,1) = this%rbuffzC(:,:,1,1)
    call transpose_z_to_y(this%rbuffzE(:,:,:,1), this%rbuffyE(:,:,:,1), this%gpE)
    call transpose_y_to_x(this%rbuffyE(:,:,:,1), this%tauijWM(:,:,:,2), this%gpE)

    this%ustarsqvar = this%WallMFactorvar * this%filteredSpeedSq(:,:,1)

    ! NOTE:: tauijWMhat_inY and tauijWMhat_inZ are not populated. Are they
    ! required ????

end subroutine

subroutine getSpanAvgVelAtWall(this, uhatC, vhatC)
    class(sgs_igrid), intent(inout), target :: this
    complex(rkind), dimension(this%sp_gpC%ysz(1),this%sp_gpC%ysz(2),this%sp_gpC%ysz(3)), intent(in) :: uhatC, vhatC

    real(rkind), dimension(:,:,:), pointer :: rbuffx1, rbuffx2, rbuffx3
    complex(rkind), dimension(:,:,:), pointer :: cbuffy, tauWallH

    cbuffy => this%cbuffyC(:,:,:,1); tauWallH => this%cbuffzC(:,:,:,1)     

end subroutine

subroutine getfilteredSpeedSqAtWall(this, uhatC, vhatC)
    class(sgs_igrid), intent(inout), target :: this
    complex(rkind), dimension(this%sp_gpC%ysz(1),this%sp_gpC%ysz(2),this%sp_gpC%ysz(3)), intent(in) :: uhatC, vhatC

    real(rkind), dimension(:,:,:), pointer :: rbuffx1, rbuffx2, rbuffx3
    complex(rkind), dimension(:,:,:), pointer :: cbuffy, tauWallH

    cbuffy => this%cbuffyC(:,:,:,1); tauWallH => this%cbuffzC(:,:,:,1)     
    rbuffx3 => this%filteredSpeedSq; rbuffx1 => this%rbuffxC(:,:,:,1)
    rbuffx2 => this%rbuffxC(:,:,:,2)

    call transpose_y_to_z(uhatC,tauWallH,this%sp_gpC)
    call this%spectC%SurfaceFilter_ip(tauWallH(:,:,1))
    call transpose_z_to_y(tauWallH,cbuffy, this%sp_gpC)
    call this%spectC%ifft(cbuffy,rbuffx1)

    call transpose_y_to_z(vhatC,tauWallH,this%sp_gpC)
    call this%spectC%SurfaceFilter_ip(tauWallH(:,:,1))
    call transpose_z_to_y(tauWallH,cbuffy, this%sp_gpC)
    call this%spectC%ifft(cbuffy,rbuffx2)

    if(this%is_z0_varying) then
        this%Uxvar = rbuffx1(:,:,1)
        this%Uyvar = rbuffx2(:,:,1)
    endif

    rbuffx1 = rbuffx1*rbuffx1
    rbuffx2 = rbuffx2*rbuffx2
    rbuffx3 = rbuffx1 + rbuffx2

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

          !if(nrank==0) then
          !   write(nrank+100,'(8(e19.12,1x),2(i5,1x))') this%ustar, this%invObLength, this%Tsurf, this%wTh_surf, ustarDiff, this%PsiM, u, PsiH, idx, itermax
          !endif
          ! Inside the do loop all the used variables are on the stored on the stack
          ! After the while loop these variables are copied to their counterparts
          ! on the heap (variables part of the derived type)
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
          !if(nrank==0) then
          !   write(nrank+200,'(8(e19.12,1x),2(i5,1x))') this%ustar, this%invObLength, this%Tsurf, this%wTh_surf, ustarDiff, this%PsiM, u, PsiH, idx, itermax
          !endif
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
   
          ! Inside the do loop all the used variables are on the stored on the stack
          ! After the while loop these variables are copied to their counterparts
          ! on the heap (variables part of the derived type)
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
      if(this%is_z0_varying) then
          ! computed in BouZeidLocalModel
      else
          this%ustar = this%Uspmn*kappa/(log(hwm/this%z0))
      endif
      this%invObLength = zero
      this%wTh_surf = zero
      this%PsiM = zero
    end if 
end subroutine
