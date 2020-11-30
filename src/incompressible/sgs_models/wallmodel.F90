subroutine destroyWallModel(this)
   class(sgs_igrid), intent(inout) :: this
   deallocate(this%tauijWM, this%tauijWMhat_inZ, this%tauijWMhat_inY)
   if (allocated(this%filteredSpeedSq)) deallocate(this%filteredSpeedSq)
   !if (allocated(this%mask_upstream)) deallocate(this%mask_upstream)
end subroutine

subroutine initWallModel(this)
   class(sgs_igrid), intent(inout) :: this
   real(rkind) :: epssmall = 1.0d-6

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
      if(this%is_z0_varying) then
          this%kaplnzfac_s = (kappa/log(half*this%dz/this%z0s))**2
          this%kaplnzfac_r = (kappa/log(half*this%dz/this%z0r))**2
      endif
   case (3) ! Abkar-PA heterogeneous Wall model
      if(.not. this%is_z0_varying) then
          call GracefulExit("You cannot use Abkar-PA wall model for a homogeneous problem", 111)
      endif
      this%kaplnzfac_s = (kappa/log(half*this%dz/this%z0s))**2
      !if(nrank==0) then
      !   write(*,'(a,5(1x,e19.12))') 'kaplnzfac_s:= ', kappa, this%dz, this%z0s, this%kaplnzfac_s
      !endif
      this%kaplnzfac_r = (kappa/log(half*this%dz/this%z0r))**2
      !allocate(this%mask_upstream(this%gpC%xsz(1), this%gpC%xsz(2)))
      allocate(this%filteredSpeedSq(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)))
      !where(this%lamfact > (one-epssmall))
      !    this%mask_upstream = one
      !elsewhere
      !    this%mask_upstream = zero
      !endwhere
      !this%mask_normfac = p_sum(sum(this%mask_upstream))
   case (4) ! Our heterogeneous Wall model
      this%kaplnzfac_s = (kappa/log(half*this%dz/this%z0s))**2
      this%kaplnzfac_r = (kappa/log(half*this%dz/this%z0r))**2
      !allocate(this%mask_upstream(this%gpC%xsz(1), this%gpC%xsz(2)))
      allocate(this%filteredSpeedSq(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)))
      ! condition mask_upstream on deli and not on lamfact
      !where(this%deli < epssmall)
      !    this%mask_upstream = one
      !elsewhere
      !    this%mask_upstream = one
      !endwhere
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
   real(rkind) :: ust1fac, ustar1, epssmall = 1.0d-6, dzby2
   integer, dimension(this%gpC%xsz(1), this%gpC%xsz(2)) :: modelregion
   integer :: i, j

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
      !! type 1 :: SG-local method (span averages)
      !call this%getSpanAvgVelAtWall(uhat, vhat)
      !call this%compute_ustar_upstreampart(ustar1)

      ! type 2 :: BZ-local method (twice-filtered)
      if(this%filter_for_heterog) then
          call this%getfilteredSpeedSqAtWall(uhat, vhat)
      else
          if(this%gpC%xst(3)==1) then
              this%Uxvar = u(:,:,1)
              this%Uyvar = v(:,:,1)
              this%filteredSpeedSq = u*u + v*v
          endif
      endif
      call this%getSpanAvgVelAtWall()
      this%rbuffxC(:,:,1,2) = sqrt(this%filteredSpeedSq(:,:,1))

      ! using this%filteredSpeedSq in the upstream region, estimate ustar1
      call this%compute_ustar_upstreampart(ustar1)
      ust1fac = ustar1/sqrt(this%kaplnzfac_s)
      this%ustar_upstream = ustar1

      !modelregion = 0; 
      where(this%lamfact > (one-epssmall))
          this%ustarsqvar = this%kaplnzfac_s*this%filteredSpeedSq(:,:,1)
          !modelregion = 1
      elsewhere (this%lamfact > epssmall)
          this%ustarsqvar = (this%rbuffxC(:,:,1,2) - this%lamfact*ust1fac) / (one - this%lamfact)
          this%ustarsqvar = this%kaplnzfac_r*this%ustarsqvar*this%ustarsqvar
          !modelregion = 2
      elsewhere
          this%ustarsqvar = this%kaplnzfac_r*this%filteredSpeedSq(:,:,1)
          !modelregion = 3
      endwhere

      !if(nrank==0) then
      !  !do i=1,this%gpC%xsz(1)
      !  !  write(*,*) this%lamfact(i,1), modelregion(i,1)
      !  !enddo
      !  i = 34
      !  write(*,'(a,e19.12,x,i4.4,1x,4(e19.12,1x))') '---i=34: ', this%lamfact(i,1), modelregion(i,1), this%rbuffxC(i,1,1,2), ust1fac, this%ustarsqvar(i,1), this%kaplnzfac_r
      !endif

      ! tau_13
      this%rbuffxC(:,:,1,1) = -this%ustarsqvar * this%Uxvar / (this%rbuffxC(:,:,1,2) + 1.0d-18)
      call this%set_tauijWM(this%rbuffxC(:,:,:,1), 1)

      ! tau_23
      this%rbuffxC(:,:,1,1) = -this%ustarsqvar * this%Uyvar / (this%rbuffxC(:,:,1,2) + 1.0d-18)
      call this%set_tauijWM(this%rbuffxC(:,:,:,1), 2)

   case (4) ! Our heterogeneous model

      ! how do we determine ustarsqvar??
      if(this%filter_for_heterog) then
          call this%getfilteredSpeedSqAtWall(uhat, vhat)
      else
          if(this%gpC%xst(3)==1) then
              this%Uxvar = u(:,:,1)
              this%Uyvar = v(:,:,1)
              this%filteredSpeedSq = u*u + v*v
          endif
      endif
      call this%getSpanAvgVelAtWall()
      this%rbuffxC(:,:,1,2) = sqrt(this%filteredSpeedSq(:,:,1))

      ! using this%filteredSpeedSq in the upstream region, estimate ustar1
      call this%compute_ustar_upstreampart(ustar1)
      !ust1fac = ustar1/sqrt(this%kaplnzfac_s)

      dzby2 = half*this%dz 
      do j = 1, this%gpC%xsz(2)
        do i = 1, this%gpC%xsz(1)
            if(this%deli(i,j) < epssmall) then
                this%ustarsqvar(i,j) = this%kaplnzfac_s*this%filteredSpeedSq(i,j,1)
            elseif(this%deli(i,j) < dzby2) then
                call this%solve_nonlinprob_1(i, j)
            else
              if(this%alpfac*this%deli(i,j) < dzby2) then
                  call this%solve_nonlinprob_2(i, j, ustar1)
              else
                  this%ustarsqvar(i,j) = this%kaplnzfac_r*this%filteredSpeedSq(i,j,1)
              endif
            endif
        enddo
      enddo

      ! tau_13
      this%rbuffxC(:,:,1,1) = -this%ustarsqvar * this%Uxvar / (this%rbuffxC(:,:,1,2) + 1.0d-18)
      call this%set_tauijWM(this%rbuffxC(:,:,:,1), 1)

      ! tau_23
      this%rbuffxC(:,:,1,1) = -this%ustarsqvar * this%Uyvar / (this%rbuffxC(:,:,1,2) + 1.0d-18)
      call this%set_tauijWM(this%rbuffxC(:,:,:,1), 2)

   end select

end subroutine

subroutine solve_nonlinprob_1(this, i, j) !, ustar2)
   class(sgs_igrid), intent(inout) :: this
   integer, intent(in) :: i, j

   real(rkind) :: deli, dele, z01, z02, lnfdez02, lnfdiz01, tol, xvar, xvar_new
   real(rkind) :: xdiff, term1, termPR, termr, terms, Rfac, ustar1_sq, onemxsq
   integer :: iter, max_iters

   deli = this%deli(i,j); z02 = this%z0r; z01 = this%z0s; dele = deli*this%alpfac
   lnfdez02 = log(dele/z02)/kappa
   lnfdiz01 = log(deli/z01)/kappa

   max_iters = 1000; tol = 1.0d-8 
   xvar = zero; 
   ! Set initial guess based on the type of transition
   if(z01 < z02) then
       xvar_new = 1.1_rkind   ! S-R Transition
   else
       xvar_new = 0.9_rkind   ! R-S Transition
   endif
   do iter = 1, max_iters
     xdiff = abs(one - xvar/xvar_new)
     if(xdiff < tol) exit
     xvar = xvar_new;

     onemxsq = (one-xvar)**2
     termr = -(one + (deli - xvar*dele)/(four*this%betfac*(deli+xvar*dele)))
     terms = half*xvar*dele/(this%betfac*(deli+xvar*dele))
     Rfac  = sqrt(termr*termr-four*terms)
     termPR = half*(one-xvar**2) + (deli - xvar*dele)/(four*this%betfac*(deli+xvar*dele))
     term1 = -onemxsq + (xvar*xvar + termr*termPR - onemxsq*terms)/Rfac
     term1 = term1 + termPR*log(terms/(one+termr+terms))

     xvar_new = (term1 + lnfdiz01)/lnfdez02
   enddo
   if(xdiff > tol) then
     call message(1, "Did not converge", 0.0_rkind)
     call message(2, "iter: ", iter  )
     call message(2, "xdiff: ", xdiff)
     call message(2, "xvar : ", xvar )
     call GracefulExit("Nonlinear solver in wall model did not converge. Check details.", 999)
   endif
   ustar1_sq = this%kaplnzfac_s*this%filteredSpeedSq(i,j,1)
   this%ustarsqvar(i,j) = xvar*xvar*ustar1_sq

end subroutine

subroutine solve_nonlinprob_2(this, i, j, ustar1)
   class(sgs_igrid), intent(inout) :: this
   integer, intent(in) :: i, j
   real(rkind), intent(in)  :: ustar1

   real(rkind) :: deli, dele, z01, z02, lnfdiz01, tol, xvar, xvar_new
   real(rkind) :: xdiff, onemxsq, term1, termPR, termr, terms, Rfac
   real(rkind) :: zbar, ufix
   integer :: iter, max_iters

   deli = this%deli(i,j); z02 = this%z0r; z01 = this%z0s; dele = deli*this%alpfac
   lnfdiz01 = log(deli/z01)/kappa

   zbar = (half*this%dz - dele)/(deli-dele)
   ufix = sqrt(this%filteredSpeedSq(i,j,1))

   max_iters = 1000; tol = 1.0d-8 
   xvar = zero; 
   ! Set initial guess based on the type of transition
   if(z01 < z02) then
       xvar_new = 1.1_rkind   ! S-R Transition
   else
       xvar_new = 0.9_rkind   ! R-S Transition
   endif
   do iter = 1, max_iters
     xdiff = abs(one - xvar/xvar_new)
     if(xdiff < tol) exit
     xvar = xvar_new;

     onemxsq = (one-xvar)**2
     termr = -(one + (deli - xvar*dele)/(four*this%betfac*(deli+xvar*dele)))
     terms = half*xvar*dele/(this%betfac*(deli+xvar*dele))
     Rfac  = sqrt(termr*termr-four*terms)
     termPR = half*(one-xvar**2) + (deli - xvar*dele)/(four*this%betfac*(deli+xvar*dele))
     term1 = (zbar-one) * onemxsq + (xvar*xvar + termr*termPR - onemxsq*terms)/Rfac
     term1 = term1 + termPR*log((zbar*zbar+zbar*termr+terms)/(one+termr+terms))

     xvar_new = -half*(deli-dele)*term1/(this%betfac*kappa*dele*(ufix/ustar1 - lnfdiz01))
   enddo
   if(xdiff > tol) then
     call message(1, "Did not converge", 0.0_rkind)
     call message(2, "iter: ", iter  )
     call message(2, "xdiff: ", xdiff)
     call message(2, "xvar : ", xvar )
     call GracefulExit("Nonlinear solver in wall model did not converge. Check details.", 999)
   endif
   this%ustarsqvar(i,j) = (xvar*ustar1)**2

end subroutine
 
subroutine compute_ustar_upstreampart(this, ustar1)
   class(sgs_igrid), intent(inout) :: this
   real(rkind), intent(out) :: ustar1
   real(rkind) :: ufiltavg

   !ufiltavg = p_sum(sum(sqrt(this%filteredSpeedSq(:,:,1)*this%mask_upstream)))/this%mask_normfac
   !ustar1 = ufiltavg*sqrt(this%kaplnzfac_s)

   ufiltavg = p_sum(sum(this%filteredSpeedSq(:,:,1)*this%mask_upstream))/this%mask_normfac
   ustar1 = sqrt(ufiltavg*this%kaplnzfac_s)
   !if(nrank==0) 
   !print '(a,3(e19.12,1x))', "ustar1:= ", ustar1, this%kaplnzfac_s, ufiltavg

end subroutine

subroutine set_tauijWM(this, rbuffx1, ind)
   class(sgs_igrid), intent(inout) :: this
   real(rkind), dimension(this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3)), intent(in) :: rbuffx1
   integer, intent(in) :: ind

   call transpose_x_to_y(rbuffx1, this%rbuffyC(:,:,:,1), this%gpC)
   call transpose_y_to_z(this%rbuffyC(:,:,:,1), this%rbuffzC(:,:,:,1), this%gpC)

   this%rbuffzE = 0.0d0;    this%rbuffzE(:,:,1,1) = this%rbuffzC(:,:,1,1)
   call transpose_z_to_y(this%rbuffzE(:,:,:,1), this%rbuffyE(:,:,:,1), this%gpE)
   call transpose_y_to_x(this%rbuffyE(:,:,:,1), this%tauijWM(:,:,:,ind), this%gpE)

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
    real(rkind) :: ustar1

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

    this%ustarsqvar = -this%WallMFactorvar * this%filteredSpeedSq(:,:,1)

    call this%compute_ustar_upstreampart(ustar1)
    this%ustar_upstream = ustar1
    !print *, 'ustar_upstream:= ',  this%ustar_upstream

    ! NOTE:: tauijWMhat_inY and tauijWMhat_inZ are not populated. Are they
    ! required ????

end subroutine

subroutine getSpanAvgVelAtWall(this)
    class(sgs_igrid), intent(inout), target :: this
    integer :: j, k

    ! copmute span-avg of Uxvar, Uyvar, filteresSpeedSqAtWall
    this%filteredSpeedSq(:,:,2) = this%Uxvar(:,:)
    this%filteredSpeedSq(:,:,3) = this%Uyvar(:,:)
    call transpose_x_to_y(this%filteredSpeedSq, this%rbuffyC(:,:,:,1), this%gpC)

    !this%rbuffyC(:,1,1,1) = sum(this%rbuffyC(:,:,1,1),2)/real(this%gpC%ysz(2), rkind)
    this%rbuffyC(:,1,2,1) = sum(this%rbuffyC(:,:,2,1),2)/real(this%gpC%ysz(2), rkind)
    this%rbuffyC(:,1,3,1) = sum(this%rbuffyC(:,:,3,1),2)/real(this%gpC%ysz(2), rkind)
    do k=2,3
      do j=2,this%gpC%ysz(2)
          this%rbuffyC(:,j,k,1) = this%rbuffyC(:,1,k,1)
      enddo
    enddo

    call transpose_y_to_x(this%rbuffyC(:,:,:,1), this%filteredSpeedSq, this%gpC)
    this%Uxvar(:,:) = this%filteredSpeedSq(:,:,2)
    this%Uyvar(:,:) = this%filteredSpeedSq(:,:,3)

    this%filteredSpeedSq(:,:,1) = this%filteredSpeedSq(:,:,2)*this%filteredSpeedSq(:,:,2)
    this%filteredSpeedSq(:,:,1) = this%filteredSpeedSq(:,:,1) + this%filteredSpeedSq(:,:,3)*this%filteredSpeedSq(:,:,3)
    
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

    if(this%is_z0_varying .and. (this%gpC%xst(3)==1)) then
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
