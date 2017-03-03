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
      this%tauijWMhat_inZ(:,:,1,1) = this%WallMFactor*cbuffz(:,:,1) 
      call transpose_z_to_y(this%tauijWMhat_inZ(:,:,:,1), this%tauijWMhat_inY(:,:,:,1), this%sp_gpE)
      call this%spectE%ifft(this%tauijWMhat_inY(:,:,:,1), this%tauijWM(:,:,:,1))

      ! Tau_23
      call transpose_y_to_z(vhat, cbuffz, this%sp_gpC)
      this%tauijWMhat_inZ(:,:,1,2) = this%WallMFactor*cbuffz(:,:,1) 
      call transpose_z_to_y(this%tauijWMhat_inZ(:,:,:,2), this%tauijWMhat_inY(:,:,:,2), this%sp_gpE)
      call this%spectE%ifft(this%tauijWMhat_inY(:,:,:,2), this%tauijWM(:,:,:,2))

   case (2) ! Bou-zeid Wall model 
      this%WallMFactor = -(kappa/(log(this%dz/(two*this%z0)) + beta_m*this%InvObLength*this%dz/two))**2 
      call this%getfilteredSpeedSqAtWall(uhat, vhat)
      
      call this%spectC%fft(this%filteredSpeedSq, cbuffy)
      call transpose_y_to_z(cbuffy, cbuffz, this%sp_gpC)
      
      ! tau_13
      this%tauijWMhat_inZ(:,:,1,1) = (this%WallMFactor*this%umn/this%Uspmn) * cbuffz(:,:,1) 
      call transpose_z_to_y(this%tauijWMhat_inZ(:,:,:,1), this%tauijWMhat_inY(:,:,:,1), this%sp_gpE)
      call this%spectE%ifft(this%tauijWMhat_inY(:,:,:,1), this%tauijWM(:,:,:,1))
      
      ! tau_23
      this%tauijWMhat_inZ(:,:,1,2) = (this%WallMFactor*this%vmn/this%Uspmn) * cbuffz(:,:,1) 
      call transpose_z_to_y(this%tauijWMhat_inZ(:,:,:,2), this%tauijWMhat_inY(:,:,:,2), this%sp_gpE)
      call this%spectE%ifft(this%tauijWMhat_inY(:,:,:,2), this%tauijWM(:,:,:,2))
   end select

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
        this%Umn = real(uhat(1,1,1),rkind)*this%meanFact
        this%Vmn = real(vhat(1,1,1),rkind)*this%meanFact
        this%Uspmn = real(cbuff(1,1,1),rkind)*this%meanFact
        if (this%isStratified) this%Tmn = real(That(1,1,1),rkind)*this%meanFact
    end if
    call mpi_bcast(this%Umn,1,mpirkind,0,mpi_comm_world,ierr)
    call mpi_bcast(this%Vmn,1,mpirkind,0,mpi_comm_world,ierr)
    call mpi_bcast(this%Uspmn,1,mpirkind,0,mpi_comm_world,ierr)
    if (this%isStratified) call mpi_bcast(this%Tmn,1,mpirkind,0,mpi_comm_world,ierr)

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

subroutine getSurfaceQuantities(this)
    class(sgs_igrid), intent(inout) :: this
    integer :: idx
    integer, parameter :: itermax = 100 
    real(rkind) :: ustarNew, ustarDiff, dTheta, ustar
    real(rkind) :: a, b, c, PsiH, PsiM, wTh, z, u, Linv
  
    select case (this%botBC_Temp)
    case(0) ! Dirichlet BC for temperature 
        dTheta = this%Tsurf - this%Tmn; Linv = zero
        z = this%dz/two ; ustarDiff = one; wTh = zero
        a=log(z/this%z0); b=beta_h*this%dz/two; c=beta_m*this%dz/two 
        PsiM = zero; PsiH = zero; idx = 0; ustar = one; u = this%Uspmn
   
        ! Inside the do loop all the used variables are on the stored on the stack
        ! After the while loop these variables are copied to their counterparts
        ! on the heap (variables part of the derived type)
        do while ( (ustarDiff > 1d-12) .and. (idx < itermax))
            ustarNew = u*kappa/(a - PsiM)
            wTh = dTheta*ustarNew*kappa/(a - PsiH) 
            Linv = -kappa*wTh/((this%Fr**2) * this%ThetaRef*ustarNew**3)
            PsiM = -c*Linv; PsiH = -b*Linv;
            ustarDiff = abs((ustarNew - ustar)/ustarNew)
            ustar = ustarNew; idx = idx + 1
        end do 
        this%ustar = ustar; this%invObLength = Linv; this%wTh_surf = wTh
    case(1) ! Homogeneous Neumann BC for temperature
        this%ustar = this%Uspmn*kappa/(log(this%dz/two/this%z0))
        this%invObLength = zero
        this%wTh_surf = zero
    end select
end subroutine
